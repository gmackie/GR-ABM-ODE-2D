Windows 7 Build Instructions {#readme-win7}
============================

Cygwin
------

-# Go to http://cygwin.com/install.html
-# Download the setup.exe file listed there and run it.
-# Follow the installation instructions until you see a list of packages
   - *Try to pick a download site near you (like ftp://lug.mtu.edu)*
-# Select the following packages by clicking on circle arrow icon to get all the required libraries/utilities
   - To get started, you will need the following packages:
     + subversion
     + gcc4-g++
     + libboost-devel
     + make
   - To compile the gui, you will need the following packages
     + xinit
     + freeglut
     + qt4-devel-tools
     + libGLU-devel
     + libQtOpenGL4-devel
     + libQtSvg4-devel
     + libQtXml4-devel
     + libQtDesigner4-devel
-# Click next, then finish

**Note: the following steps are needed only to compile the gui**
-# Open the cygwin terminal application (most likely on the desktop)
-# Type the following:

        $ svn co https://qwt.svn.sourceforge.net/svnroot/qwt/branches/qwt-6.0
        $ cd qwt-6.0
        $ qmake-qt4
        $ make && make install
        $ cp qwt.prf qwtconfig.pri /usr/share/qt4/mkspecs/features

Our code
--------

**Note:** currently the code expects the boost libraries to be of the form libboost_library_name.a
This means in Cygwin, you'll get a linker error at the end of the compilation (until we fix it in the makefile that is) due to the extra -mt at the end.
To deal with this, change all -lboost_library_name to -lboost_library_name-mt within the makefile and grviz.pro
- Open Cygwin/X application and type the following:

      $ svn checkout svn://innoculant.micro.med.umich.edu/dev/gr2d/GR-ABM-ODE path/to/my/repo
      $ cd path/to/my/repo
      $ cd simulation
      $ make
      $ ./gr --help
- For the gui version, type the following:

      $ cd path/to/my/repo
      $ qmake-qt4 -spec cygwin-g++
      $ make
      $ ./grviz-lung
