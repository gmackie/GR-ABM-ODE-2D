/*
 * Serialization.h
 *
 *  Created on: Jun 17, 2011
 *      Author: pwolberg
 */

#ifndef SERIALIZATION_H_
#define SERIALIZATION_H_

#include <iostream>
#include <string>

/*
 * Support functions for serialization/deserialization.
 *
 * The simulation code has the capability to save and load simulation
 * states.
 *
 * A loaded simulation state can be used in the GUI version of
 * the code to display and save pictures for that state. The non-GUI
 * version can be run on a cluster (which generally don't have graphics
 * support - a lot of CPUs but no graphics cards) and set to save states
 * periodically. Once the all the cluster runs (typically in the hundreds
 * or thousands) have finished a script can be used to post process the
 * saved states, using the GUI version to load each state for each run,
 * saving pictures for that state. This allows looking at both numerical
 * statistics and graphics for each run, avoiding the need to make a decision
 * about which runs are interesting on the basis of statistics alone and
 * rerunning just those runs using the GUI version on a non-cluster environment
 * to generate graphics for just those runs.
 *
 * Another use of saved states is to re-run a simulation to get more detail
 * about its results, perhaps looking at graphics or statistics at more frequent
 * intervals (every day rather than every 50 days or every time step rather then
 * every day) for some key part of the simulation. Using a saved state allows
 * re-starting a simulation near the point of interest rather than from the beginning,
 * saving considerable time, especially for those models, like 3D models, that run
 * rather slowly.
 *
 * Loading a saved state is useful only if it accurately represents the simulation state.
 * Given the same random number generator seed (from a saved random number generator
 * state) and the same parameter file, a loaded state should produce exactly the same
 * graphics and simulation results from the saved state onward as the original simulation
 * run that produced the saved state. For example, if a simulation is run for 200 simulation
 * days and its state is saved at day 100, with graphics saved every 50 days and statistics
 * saved every day. If the state saved at day 100 is loaded and rerun until day 200, with
 * pictures and statistics saved at the same time intervals, the pictures at day 100, 150
 * and 200 from the run of the loaded state should be exactly the same as from the original
 * run. Similarly the statistics for day 100 to 200 from the run of the loaded state should
 * be exactly the same as from the original run.
 *
 * The only caveat here is that the random number generators aren't guaranteed to work the same
 * on different systems, or with different versions of the random number library. But on the same
 * system or type of system (ex. Linux on a PC system), with the same version of the random number
 * library and with the same parameter file, the results should be the same.
 *
 * Also note that we are saving states in text form. This means that they are human readable and also
 * that this avoids the differences between 32 bit and 64 bit systems (ex. different sizes for ints, etc.)
 * and big endian vs. little endian issues.
 *
 * The classes in the simulation code that represent the simulation state have
 * serialization and deserialization functions. These functions must be updated to reflect any
 * changes in the members for a class, otherwise problems will occur when loading and running
 * a saved state.
 *
 * For example, suppose a class has a member added, but the serialization
 * and deserialization functions for that class are not updated. If the new
 * version of the class is serialized and then deserialized the new member
 * won't be saved and restored. It may have a default value defined
 * in a constructor or it may have some random value. This means a saved simulation state
 * when loaded and run won't reproduce the exact same results as when that
 * simulation was originally run (assuming the same seed and parameter file are used).
 * Unfortunately this particular problem is not addressed by the support functions here.
 * Doing so would probably require RTTI (run time type information), which is more
 * complexity and run time overhead than we want to incur. The only approach to dealing
 * with this problem is to take care when changing a class's members to update any
 * serialization/deserialization code in that class and to test that running a simulation
 * with a saved simulation state produces the same results as running the simulation
 * from the initial conditions, with the same random number seed and parameter file.
 *
 * Deleting a member generally isn't a problem because serialization code that writes
 * that member will no longer compile.
 *
 * Changing a member's type, for members that are of a primitive type (int, float, etc.),
 * generally isn't a problem because the stream output and input operators are used (">>" and "<<"),
 * which can handle any primitive type.
 * Changing a member's type from one class to another generally isn't a problem because we typically
 * call the class's serialize and deserialize functions. These function calls are (should be) the
 * same for all classes that we need to include in a saved state. That is, if class Flip to be serialized
 * has a member of type Class Gorp which is changed to be of type Class Blap, both Gorp and Blap
 * should have functions named "serialize" and "deserialize", so changing that member's type shouldn't
 * require any changes to the Flip class serialize and deserialize functions.
 *
 * The support functions in this class that write a header and footer at the start and end of the data
 * written for a class help protect against a mismatch between what a class's serialization function
 * writes and what its deserialization function reads. It also makes it easier for a person to determine
 * the start and end of the data for a serialized class when reading a saved state.
 *
 * The header and footer should be single words - no white space. This makes reading them easy using a
 * single input stream operation (i.e. in << text) rather than having to read several pieces of text.
 *
 * If a class member is written by a class's serialization function but not read by its deserialization
 * function then the number of items written will be larger than the number that are read.
 * The deserialization function will attempt to read the footer for the class before the footer is reached, so
 * what it reads as the footer won't match its actual footer and be reported as an error.
 *
 * If a class member is not written by a class's serialization function but is read by its deserialization
 * function then the number of items written will be smaller than the number that are read.
 * The deserialization function will read data after the footer for the class being read (for the next class
 * that was serialized). Since the footer is text, if the member that the footer is read into is not a string
 * then an error will occur and be reported. If that doesn't happen then later what it reads as the footer won't
 * match the actual footer for the class being read and be reported as an error. If the class being read is
 * the last class for a saved state it will attempt to read past end of file which does not seem to be detected
 * as an error. At least, 0 is assigned to numeric variables being read.
 *
 * Also note that if a floating point value being serialized has a value of nan (not a number) what the output
 * stream operator (">>") writes is the text "nan". This will typically not deserialize properly, causing
 * a mismatched footer error. The solution is to find what caused a nan to be defined and written and to
 * fix that problem.
 *
 */

class Serialization
{
public:
  Serialization();
  virtual ~Serialization();

  static const std::string _HeaderSuffix;
  static const std::string _FooterSuffix;

  static void writeHeader(std::ostream& out, std::string className);
  static void writeFooter(std::ostream& out, std::string className);
  static bool readHeader(std::istream& in, std::string className);
  static bool readFooter(std::istream& in, std::string className);


protected:

  static bool readHeaderFooter(std::istream& in, std::string className, std::string header, std::string type);

};

#endif /* SERIALIZATION_H_ */
