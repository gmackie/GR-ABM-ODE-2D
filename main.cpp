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

unsigned long getSeed()
{
	const QTime curTime = QTime::currentTime();
	return curTime.hour() * 3600 + curTime.minute() * 60 + curTime.second();
}

void setupScriptingMode(MainWindow& mainWindow)
{
	Ui::MainWindowClass& ui = mainWindow.getUI();

	ui.pushButtonAnimation->click();
	ui.checkBoxDrawSmoke->setChecked(false);
	ui.groupBoxOutput->setEnabled(false);
	ui.groupBoxStoppingCriteria->setEnabled(false);
	ui.checkBoxDrawGranulomaBorder->setChecked(true);
	ui.horizontalSliderGranulomaBorderThreshold->setEnabled(false);

	mainWindow.hide();
}

Snapshot* setupOutput(MainWindow& mainWindow, const std::string& outputDir)
{
	Ui::MainWindowClass& ui = mainWindow.getUI();
	QString dirName = QString(outputDir.c_str()) + QDir::separator();
	QString fileName = QString("%1%2seed%3.csv").arg(dirName).arg(QDir::separator()).arg(g_Rand.getSeed());
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
	unsigned long seed;
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
	int nDays;
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
		("seed,s", po::value<unsigned long>(&seed))
		("diffusion,d", po::value<int>(&diffMethod)->default_value(0),
				"Diffusion method:\n0 - FTCS\n1 - BTCS (SOR, correct)\n2 - BTCS (SOR, wrong)")
		("days", po::value<int>(&nDays)->default_value(200), "Number of days to simulate")
		("script,c", "Scripting mode")
		("output,o", po::value<std::string>(&outputDir), "Output directory")
		("csv-interval", po::value<int>(&csvInterval)->default_value(1),
				"CSV update interval (10 min timesteps)")
		("png-interval", po::value<int>(&pngInterval)->default_value(144*50),
				"PNG snapshot interval (10 min timesteps)")
		("state-interval", po::value<int>(&stateInterval)->default_value(144*50),
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
		("load-state,l",  po::value<std::string>(&stateFileName), "File name of saved state to load")
		("ode", "Use integrated lymph node ODE for recruitment")
		("ln-ode", po::value<std::string>(&lymphNodeODE), "Lymph node application")
		("ln-ode-temp", po::value<std::string>(&lymphNodeTemp), "Lymph node temp file")
		("version,v", "Version number");

	try
	{
		po::positional_options_description p;
		p.add("input-file", -1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);

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
			seed = getSeed();
		}
		g_Rand.setSeed(seed);

		scriptingMode = vm.count("script");
		outputEnabled = vm.count("output");

		ode = vm.count("ode");

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

		if (!Params::getInstance(ode)->fromXml(inputFileName.c_str()))
			return 1;

		if (!(0 <= diffMethod && diffMethod < 3))
		{
			printUsage(argv[0], desc);
			return 1;
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


	ScalarAgentGrid agentGrid;
	AgentsVisualization agentsVisualization(Simulation::_DIM, &agentGrid);
	MainInterface itfc(&agentsVisualization, &agentGrid);
	GLWindow glWindow(&itfc);
	ParamWindow paramWindow(&itfc);
	MainWindow w(&itfc, &glWindow, &paramWindow, new StatWidget(), new AgentsWidget(&agentsVisualization));

	/* set recruitment method */
	if (ode)
		itfc.getSimulation().setRecruitment(new RecruitmentLnODEPure());
	else if (lymphNodeODE == "" && lymphNodeTemp == "")
		itfc.getSimulation().setRecruitment(new RecruitmentProb());
	else
		itfc.getSimulation().setRecruitment(new RecruitmentLnODE(lymphNodeODE, lymphNodeTemp));


	glWindow.resizeGLWidget(resWidth, resHeight);
	Ui::MainWindowClass& ui = w.getUI();

	ui.comboBoxDiffusion->setCurrentIndex(diffMethod);
	ui.spinBoxStopDays->setValue(nDays);

	ui.spinBoxSnapshotCsvInterval->setValue(csvInterval);
	ui.spinBoxSnapshotPicInterval->setValue(pngInterval);
	ui.spinBoxSnapshotStateInterval->setValue(stateInterval);

	if (granvizEnabled)
	{
		ui.checkBoxDrawGranulomaBorder->setChecked(granvizEnabled);
		ui.comboBoxGranulomaDataset->setCurrentIndex(granvizDatasetIndex);
	}

	itfc.getSimulation().setOutcomeMethod(0, outcomeMethod[0],
		_OUTCOME_ALPHA, testPeriod[0], samplePeriod[0]);
	itfc.getSimulation().setOutcomeMethod(1, outcomeMethod[1],
		_OUTCOME_ALPHA, testPeriod[1], samplePeriod[1]);
	w.updateOutcomeParameters();

	// Since this loads a saved seed this must come before the seed is used,
	// for example when setting up output for saving statistics.
	if (stateFileName.size() > 0)
	{
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
		pSnapshot = setupOutput(w, outputDir);
	}

	if (scriptingMode)
	{
		w.setScriptingMode(true);
		glWindow.show();
		setupScriptingMode(w);

		int err = a.exec();

		if (outputEnabled)
		{
			// take final snapshot when in scripting mode
			int time = itfc.getSimulation().getTime();
			const GrStat& stats = itfc.getSimulation().getStats();

			if (pngInterval != 0)
				pSnapshot->takePicture(time, glWindow.grabFrameBuffer());

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
		glWindow.show();
		int err = a.exec();
		return err;
	}
}
