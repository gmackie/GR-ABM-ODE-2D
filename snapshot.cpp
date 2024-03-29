/*
 * snapshot.cpp
 *
 *  Created on: Dec 3, 2009
 *      Author: mohammed
 */

#include "snapshot.h"
#include "simulation.h"
#include <QPainter>
#include <QPen>
#include <QFont>
#include <QImage>
#include <QDir>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

Snapshot::Snapshot(const QString& dirName, const QString& fileName)
  : _outFile(fileName.toLatin1().data(), std::ios_base::trunc)
  , _dirName(dirName)
{
  assert(_outFile.good());
  _outFile << "\"time\""
           << ','
           << "\"Mac\""
           << ','
           << "\"Mr\""
           << ','
           << "\"Mi\""
           << ','
           << "\"Mci\""
           << ','
           << "\"Ma\""
           << ','
           << "\"Md\""
           << ','
           << "\"Tgam\""
           << ','
           << "\"Tgam a\""
           << ','
           << "\"Tgam reg\""
           << ','
           << "\"Tgam d\""
           << ','
           << "\"Tcyt\""
           << ','
           << "\"Tcyt a\""
           << ','
           << "\"Tcyt reg\""
           << ','
           << "\"Tcyt d\""
           << ','
           << "\"Treg\""
           << ','
           << "\"Treg r\""
           << ','
           << "\"Treg d\""
           << ','
           << "\"Int. Mtb.\""
           << ','
           << "\"Ext. Mtb.\""
           << ','
           << "\"NonRepl Ext. Mtb.\""
           << ','
           << "\"Tot Mtb.\""
           << ','
           << "\"TNF\""
           << ','
           << "\"CCL2\""
           << ','
           << "\"CCL5\""
           << ','
           << "\"CXCL9\""
           << ','
           << "\"AreaTNF\""
           << ','
           << "\"AreaCellDensity\""
           << ','
           << "\"MDC\""
           << ','
           << "\"N4\""
           << ','
           << "\"TH0\""
           << ','
           << "\"TH1\""
           << ','
           << "\"N8\""
           << ','
           << "\"T80\""
           << ','
           << "\"T8\""
           << ','
           << "\"TC\""
           << ','
           << "\"TH0lung\""
           << ','
           << "\"TH1lung\""
           << ','
           << "\"T80lung\""
           << ','
           << "\"T8lung\""
           << ','
           << "\"TClung\""
           << ','
           << "\"NrSourcesMac\""
           << ','
           << "\"NrSourcesTgam\""
           << ','
           << "\"NrSourcesTcyt\""
           << ','
           << "\"NrSourcesTreg\""
           << ','
           << "\"NrCaseated\""
           << ','
           << "\"MacApoptosisTNF\""
           << ','
           << "\"MrApoptTNF\""
           << ','
           << "\"MiApoptTNF\""
           << ','
           << "\"MciApoptTNF\""
           << ','
           << "\"MaApoptTNF\""
           << ','
           << "\"TcellApoptTNF\""
           << ','
           << "\"MrActivationTNF\""
           << ','
           << "\"MiActivationTNF\""
           << "\n";
}

Snapshot::~Snapshot()
{
  if (_outFile)
    {
      _outFile.flush();
      _outFile.close();
    }
}

void Snapshot::takeStateSnapshot(int time, const Simulation& sim)
{
  int days, hours, minutes;
  GrSimulation::convertSimTime(time, days, hours, minutes);

  QString fileName = _dirName + QDir::separator() +
                     QString("%1d%2h%3m.state.gz").arg(days, 3, 10, QChar('0')).
                     arg(hours, 2, 10, QChar('0')).arg(minutes, 2, 10, QChar('0'));

  namespace bio = boost::iostreams;
  bio::filtering_ostream out;
  out.push(bio::gzip_compressor());
  out.push(bio::file_sink(fileName.toStdString()));

  if (out.good())
    {
      sim.saveState(out);
    }
  else
    {
      std::cerr << " Snapshot::takeStateSnapshot, unable to save state to " << fileName.toStdString() << std::endl;
    }
}

void Snapshot::takePicture(const int time, const QImage& image, int slice, const QString prefix)
{
  int days, hours, minutes;
  GrSimulation::convertSimTime(time, days, hours, minutes);

  QString fileName = _dirName + QDir::separator() + prefix +
                     QString("%2d%3h%4m").arg(days, 3, 10, QChar('0')).
                     arg(hours, 2, 10, QChar('0')).arg(minutes, 2, 10, QChar('0'))
                     + (slice > -1 ? QString("-slice%4.png").arg(slice, 3, 10, QChar('0')) : ".png");

  image.save(fileName);
}
void Snapshot::takePicture(const int time, const QImage& image, const QString prefix)
{
  takePicture(time, image, -1, prefix);
}

void Snapshot::takeSnapshot(const int time, const Stats& stats)
{
  _outFile << time;

  _outFile << (stats.getNrOfMacs()) << ',';
  for(size_t i=0; i<Mac::NSTATES; i++)
    _outFile << (stats.getNrOfMacs(Mac::State(i))) << ',';
  _outFile << (stats.getNrOfTgams()) << ',';
  for(size_t i=0; i<Tgam::NSTATES; i++)
    _outFile << (stats.getNrOfTgams(Tgam::State(i))) << ',';
  _outFile << (stats.getNrOfTcyts()) << ',';
  for(size_t i=0; i<Tcyt::NSTATES; i++)
    _outFile << (stats.getNrOfTcyts(Tcyt::State(i))) << ',';
  _outFile << (stats.getNrOfTregs()) << ',';
  for(size_t i=0; i<Treg::NSTATES; i++)
    _outFile << (stats.getNrOfTregs(Treg::State(i))) << ',';
  _outFile
      << stats.getTotIntMtb() << ','
      << stats.getTotExtMtb() << ','
      << stats.getTotNonRepExtMtb() << ','
      << (stats.getTotIntMtb() + stats.getTotExtMtb()) << ','
      << stats.getTotTNF() << ','
      << stats.getTotCCL2() << ','
      << stats.getTotCCL5() << ','
      << stats.getTotCXCL9() << ','
      << stats.getAreaTNF() << ','
      << stats.getAreaCellDensity() << ','
      << stats.getMDC() << ','
      << stats.getN4() << ','
      << stats.getTH0() << ','
      << stats.getTH1()	<< ','
      << stats.getN8()	<< ','
      << stats.getT80()	<< ','
      << stats.getT8()	<< ','
      << stats.getTC()	<< ','
      << stats.getTH0lung()	<< ','
      << stats.getTH1lung()	<< ','
      << stats.getT80lung()	<< ','
      << stats.getT8lung()	<< ','
      << stats.getTClung();

  int startState = 1; // Skip the dead state: apoptosis doesn't occur for an already dead mac.
  static int totMacApoptosisTNF[Mac::NSTATES] = {0}; //Just temporary for Mohammed Fallahi
  int sumMacApoptosisTNF = 0;
  for(int i=startState; i<Mac::NSTATES; i++)   //Keep a running sum of deaths
    {
      totMacApoptosisTNF[i]+=(stats.getMacApoptosisTNF((Mac::State)i));
      sumMacApoptosisTNF+=totMacApoptosisTNF[i];
    }

  _outFile << (stats.getTotNrSources()) << ',';
  for(size_t i=0; i<NAGENTS; i++)
    _outFile << (stats.getNrSources((AgentType)i)) << ',';

  _outFile << stats.getNrCaseated()
           << ','
           << sumMacApoptosisTNF
           << ',';

  for(int i=startState; i<Mac::NSTATES; i++)
    _outFile<<totMacApoptosisTNF[i]<<',';

  _outFile
      << stats.getTcellApoptosisTNF();


  _outFile << std::endl;

}
