#include "gui/mainwindow.h"
#include "gui/glwindow.h"
#include "gui/paramwindow.h"
#include "gui/statwidget.h"
#include "gui/agentswidget.h"
#include "visualization/agentsvisualization.h"
#include "scalardatasets/scalaragentgrid.h"

#include <QtGui>
#include <QApplication>
#include <QTime>
#include <QDir>
#include <QString>
#include <iostream>
#include <boost/program_options.hpp>
#include <string>
#include "simulation/params.h"
#include "simulation/recruitmentprob.h"
#include "simulation/recruitmentlnode.h"
#include "simulation/recruitmentlnodepure.h"
#include "simulation/recruitmentlnodeproxy.h"
#include "maininterface.h"
#include "simulation.h"
#include "ui_mainwindow.h"
#include "snapshot.h"

namespace po = boost::program_options;

void printUsage(char* pArgv0, po::options_description& desc)
{
	std::cout << "Usage: " << pArgv0 << " [options]\n" << desc << std::endl;
}

void printVersion()
{
	std::cout << "Version: " << GR_VERSION << std::endl;
}

void setupScriptingMode(MainWindow& mainWindow)
{
	Ui::MainWindowClass& ui = mainWindow.getUI();

	ui.pushButtonAnimation->click();
	ui.checkBoxDrawSmoke->setChecked(false);
	ui.groupBoxOutput->setEnabled(false);
	ui.groupBoxStoppingCriteria->setEnabled(false);
	ui.horizontalSliderGranulomaBorderThreshold->setEnabled(false);

	mainWindow.hide();
}

Snapshot* setupOutput(MainWindow& mainWindow, const std::string& outputDir, bool outputcsv)
{
	Ui::MainWindowClass& ui = mainWindow.getUI();
	QString dirName = QString(outputDir.c_str()) + QDir::separator();
	QString fileName = outputcsv ? QString("%1%2seed%3.csv").arg(dirName).arg(QDir::separator()).arg(g_Rand.getSeed()) : QString();
	Snapshot* pSnapshot = new Snapshot(dirName, fileName);

	mainWindow.setSnapshot(pSnapshot);

	ui.lineEditOutput->setText(dirName);
	ui.checkBoxOutput->setChecked(true);

	return pSnapshot;
}

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	bool ode;
	bool tnfrDynamics;
	bool nfkbDynamics;
	int tnfDepletionTimeStep;
	int il10DepletionTimeStep;
	unsigned int seed;
  size_t dim;
	bool seedSpecified;
	int diffMethod;
	std::string inputFileName;
	std::string outputDir;
	std::string resolution;
	std::string outcomeTest[NOUTCOMES];
	std::string lymphNodeODE;
	std::string lymphNodeTemp;
	bool scriptingMode;
	bool outputEnabled;
	int nDays, timesteps;
	int pngInterval;
	int csvInterval;
	int stateInterval;
	int resWidth, resHeight;
	int samplePeriod[NOUTCOMES] = { _OUTCOME_SAMPLE_PERIOD, _OUTCOME_SAMPLE_PERIOD };
	int testPeriod[NOUTCOMES] = { _OUTCOME_TEST_PERIOD, _OUTCOME_TEST_PERIOD };
	OutcomeMethod outcomeMethod[NOUTCOMES] = { OUTCOME_NONE, OUTCOME_NONE };
	const char* testName[NOUTCOMES] = { "test1", "test2" };

	// For enabling granuloma visualization and setting the dataset to use from the command line.
	std::string granvizDataSetName;
	bool granvizEnabled = false;
	int granvizDatasetIndex = 0;
	std::string dataSetNames;
	MainWindow::getScalarDataSetNames(dataSetNames);
	std::string argHelp = "granuloma visualization: \n" + dataSetNames;

	// For loading a saved state.
	std::string stateFileName;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "Help message")
		("input-file,i", po::value<std::string>(&inputFileName), "Input file name")
		("seed,s", po::value(&seed))
    ("dim,d", po::value(&dim)->default_value(100))
		("diffusion", po::value<int>(&diffMethod)->default_value(3),
     "Diffusion method:\n0 - FTCS\n1 - BTCS (SOR, correct)\n2 - BTCS (SOR, wrong)\n3 - FTCS Grid Swap\n4 - ADE Grid Swap")
		("timesteps,t", po::value<int>(&timesteps), "Number of time steps to simulate\nTakes precedence over --days")
		("days", po::value<int>(&nDays)->default_value(200), "Number of days to simulate")
		("script,c", "Scripting mode")
		("output,o", po::value<std::string>(&outputDir), "Output directory")
		("csv-interval", po::value<int>(&csvInterval)->default_value(0),
				"CSV update interval (10 min timesteps)")
		("png-interval", po::value<int>(&pngInterval)->default_value(0),
				"PNG snapshot interval (10 min timesteps)")
		("state-interval", po::value<int>(&stateInterval)->default_value(0),
				"State snapshot interval (10 min timesteps)")
		(testName[0], po::value<std::string>(&outcomeTest[0]),
				"Enable outcome testing based on the granuloma area (%1,%2,%3)\n"
				"%1: 'm' for Mtb and 'a' for area\n"
				"%2: sample period\n"
				"%3: test period")
		(testName[1], po::value<std::string>(&outcomeTest[1]),
				"Enable outcome testing based on the granuloma area (%1,%2,%3)\n"
				"%1: 'm' for Mtb and 'a' for area\n"
				"%2: sample period\n"
				"%3: test period")
		("res,r", po::value<std::string>(&resolution)->default_value("512x512"),
				"Resolution of the OpenGL window")
		("granuloma-visualization,g", po::value<std::string>(&granvizDataSetName), argHelp.c_str())
		("border", "Draw granuloma border")
		("load-state,l",  po::value<std::string>(&stateFileName), "File name of saved state to load")
		("snapshot",  po::value<std::string>(&stateFileName), "Load a saved state, save a graphic snapshot and quit.\nArgument is the saved state to load.")
		("recr", po::value<unsigned>()->default_value(0), "recruitment:\n0 - probability\n1 - lymph node ode proxy\n2 - lymph node ode pure")
		("tnfr-dynamics", "Use molecular level TNF/TNFR dynamics in the model")
		("NFkB-dynamics", "Use molecular level intracellular NFkB dynamics in the model")
		("tnf-depletion", po::value<int>(&tnfDepletionTimeStep)->default_value(-1), "The time step at which to stop secreting tnf, including by tnfr dynamics. -1: no depletion")
		("il10-depletion", po::value<int>(&il10DepletionTimeStep)->default_value(-1), "The time step at which to stop secreting il10, including by il10r dynamics. -1: no depletion")
		("ln-ode", po::value<std::string>(&lymphNodeODE), "Lymph node application")
		("ln-ode-temp", po::value<std::string>(&lymphNodeTemp), "Lymph node temp file")
		("quiet,q", "Don't show any visualization windows")
		("version,v", "Version number");

	po::variables_map vm;
	try
	{
		po::positional_options_description p;
		p.add("input-file", -1);
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

		std::cout << "GRID SIZE: NROWS: " << dim << " NCOLS: " << dim << std::endl;
		if (vm.count("version"))
		{
			printVersion();
			return 0;
		}

		if (vm.count("help"))
		{
			printUsage(argv[0], desc);
			return 0;
		}

		seedSpecified = vm.count("seed");
		if (!seedSpecified)
		{
			seed = createTimeSeed();
		}
		g_Rand.setSeed(seed);

		scriptingMode = vm.count("script");
		outputEnabled = vm.count("output");

		ode = vm.count("ode");
		tnfrDynamics = vm.count("tnfr-dynamics");
		nfkbDynamics = vm.count("NFkB-dynamics");

		if (!vm.count("input-file"))
		{
			QString fileName = QFileDialog::getOpenFileName(NULL, "Load parameters", "", "*.xml");
			if (fileName == QString::null)
			{
				std::cerr << "Input file missing" << std::endl;
				return 1;
			}
			else
			{
				inputFileName = fileName.toLatin1().data();
			}
		}

		// Must be done before making GrSimulation.
		// Also must be done before creating a lymph ODE recruitment object,
		// since the base lymph ODE class, RecruitmentLnODE, uses parameters in its constructor.
		if (!Params::getInstance(Pos(dim, dim))->fromXml(inputFileName.c_str()))
			return 1;

    switch (vm["diffusion"].as<int>())
    {
      case 0:
        diffMethod = DIFF_REC_EQ;
        break;
      case 3:
        diffMethod = DIFF_REC_EQ_SWAP;
        break;
      case 4:
        diffMethod = DIFF_ADE_SWAP;
        break;
      default:
        std::cerr<<"Unsupported Diffusion method"<<std::endl;
        printUsage(argv[0], desc);
        exit(1);
    }

		if (csvInterval < 0 || pngInterval < 0)
		{
			printUsage(argv[0], desc);
			return 1;
		}

		if (outputEnabled && !QDir(QString(outputDir.c_str())).exists())
		{
			std::cerr << "Specified output directory '" << outputDir
				<< "' does not exist" << std::endl;
			return 1;
		}

		if (sscanf(resolution.c_str(), "%dx%d", &resWidth, &resHeight) != 2)
		{
			std::cerr << "Invalid resolution specified" << std::endl;
			return 1;
		}

		if (!vm.count("timesteps"))
			timesteps = TIME_STEPS_PER_DAY * nDays;


		for (int i = 0; i < NOUTCOMES; i++)
		{
			if (vm.count(testName[i]))
			{
				char outcomeMethodChar;

				if (sscanf(outcomeTest[i].c_str(), "%c,%d,%d",
						&outcomeMethodChar, &samplePeriod[i], &testPeriod[i]) != 3)
				{
					std::cerr << "Please specify the outcome method, sample and test period (%c,%d,%d)"
						<< std::endl;
					return 1;
				}
				else if (samplePeriod <= 0 || testPeriod <= 0)
				{
					std::cerr << "Sample and test period should be positive integers"
						<< std::endl;
					return 1;
				}
				else
				{
					switch (outcomeMethodChar)
					{
					case 'a':
						outcomeMethod[i] = OUTCOME_AREA;
						break;
					case 'm':
						outcomeMethod[i] = OUTCOME_MTB;
						break;
					default:
						std::cerr << "Outcome method must be either area ('a') or mtb ('m')"
							<< std::endl;
						return 1;
					}
				}
			}
			else
			{
				outcomeMethod[i] = OUTCOME_NONE;
			}
		}

		granvizEnabled = vm.count("granuloma-visualization");
		if (granvizEnabled)
		{
			granvizDatasetIndex = MainWindow::getScalarDataSetIndex(granvizDataSetName);
			if (granvizDatasetIndex < 0)
			{
				std::cerr << "Invalid granuloma dataset name. It must be one of:" << std::endl;
				std::cerr << dataSetNames;
				return 1;
			}
		}

	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}

	ScalarAgentGrid agentGrid(dim);
	AgentsVisualization agentsVisualization(dim, &agentGrid);
	MainInterface itfc(Pos(dim, dim), &agentsVisualization, &agentGrid);
	GLWindow glWindow(&itfc);
	ParamWindow paramWindow(&itfc);
	MainWindow w(&itfc, &glWindow, &paramWindow, new StatWidget(), new AgentsWidget(&agentsVisualization));

	/* set recruitment method */
	// Parameters must be loaded, since since the base lymph ODE class, RecruitmentLnODE, uses parameters in its constructor.
	switch (vm["recr"].as<unsigned>())
	{
	case 0:
		itfc.getSimulation().setRecruitment(new RecruitmentProb());
	  break;
	case 1:
	  itfc.getSimulation().setRecruitment(new RecruitmentLnODEProxy());
	  break;
	case 2:
		itfc.getSimulation().setRecruitment(new RecruitmentLnODEPure());

	  break;
	default:
	  std::cerr << std::endl <<"Unsupported recruitment method" << std::endl << std::endl;
	  printUsage(argv[0], desc);
	  exit(1);
	}
	
	/* set TNF/TNFR and NFkB dynamics */
	/* When NFkB is turned on, tnfr dynamics will be turned on autamatically */
	itfc.getSimulation().setTnfrDynamics(tnfrDynamics || nfkbDynamics);
	itfc.getSimulation().setNfkbDynamics(nfkbDynamics);
	itfc.getSimulation().setTnfDepletionTimeStep(tnfDepletionTimeStep);
	itfc.getSimulation().setIl10DepletionTimeStep(il10DepletionTimeStep);


	glWindow.resizeGLWidget(resWidth, resHeight);
	Ui::MainWindowClass& ui = w.getUI();

	ui.comboBoxDiffusion->setCurrentIndex(diffMethod);
	
	 // Force a value changed signal to be emitted.
	 // If timesteps has a value that happens to be the same as the current value in the spinbox
	 // then no signal is emitted and the Simulation class doesn't get its simulation limit updated.
	 // It will remain whatever was set in the simulation constructor.
	ui.spinBoxStopTime->setValue(timesteps+1); 
	ui.spinBoxStopTime->setValue(timesteps);

	ui.spinBoxSnapshotCsvInterval->setValue(csvInterval);
	ui.spinBoxSnapshotPicInterval->setValue(pngInterval);
	ui.spinBoxSnapshotStateInterval->setValue(stateInterval);

	if (granvizEnabled)
	{
		ui.checkBoxDrawGranulomaBorder->setChecked(granvizEnabled);
		ui.comboBoxGranulomaDataset->setCurrentIndex(granvizDatasetIndex);
	}

	ui.checkBoxDrawGranulomaBorder->setChecked(vm.count("border"));

	itfc.getSimulation().setOutcomeMethod(0, outcomeMethod[0],
		_OUTCOME_ALPHA, testPeriod[0], samplePeriod[0]);
	itfc.getSimulation().setOutcomeMethod(1, outcomeMethod[1],
		_OUTCOME_ALPHA, testPeriod[1], samplePeriod[1]);
	w.updateOutcomeParameters();

	// Since this loads a saved seed this must come before the seed is used,
	// for example when setting up output for saving statistics.
	if (stateFileName.size() > 0)
	{
		QString qStateFileName = stateFileName.c_str();
	    if (!QFile::exists(qStateFileName))
	    {
	    	std::cerr << "Saved state " << stateFileName << " does not exist." << std::endl;
	    	exit(1);
	    }

		if (seedSpecified)
		{
			std::cerr << "Loading a saved state, seed " << seed << " ignored."<< std::endl;
		}
		w.loadState(stateFileName);
	}

	// This must come after any load of a saved state so that the loaded
	// seed will be used as part of the file name, not a random seed or
	// seed specified on the command line.
	Snapshot* pSnapshot = NULL;
	if (outputEnabled)
	{
		bool outputcsv = (csvInterval != 0);
		pSnapshot = setupOutput(w, outputDir, outputcsv);
	}

	if (scriptingMode)
	{
		w.setScriptingMode(true);
		if(!vm.count("quiet")) glWindow.show();

		setupScriptingMode(w);

		int err = a.exec();

		if (outputEnabled)
		{
			// take final snapshot when in scripting mode
			int time = itfc.getSimulation().getTime();
			const Stats& stats = itfc.getSimulation().getStats();

			if (pngInterval != 0){
				glWindow.updateWindow();  //Update the graphics window before taking the picture
				pSnapshot->takePicture(time, glWindow.grabFrameBuffer());
			}

			if (csvInterval != 0)
				pSnapshot->takeSnapshot(time, stats);

			if (stateInterval != 0)
				pSnapshot->takeStateSnapshot(time, itfc.getSimulation());
		}

		return err;
	}
	else
	{
		w.show();
		if(!vm.count("quiet")) glWindow.show();
		int err = a.exec();
		return err;
	}
}
