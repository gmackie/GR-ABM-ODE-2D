Linux Build Instructions      {#readme-linux}
========================

Ubuntu 12.04
------------

This also works for Linux Mint 12

    $ sudo apt-get install subversion g++ libboost-program-options-dev libboost-iostreams-dev libboost-serialization-dev

For the gui:

    $ sudo apt-get install libqt4-opengl-dev qt4-dev-tools libqwt-dev

Ubuntu 10.04
------------

-# Run the following:
    $ sudo apt-get install subversion g++ zlib1g-dev
-# Go to http://boost.org/ and download the latest version of boost_xx.bz2
-# Run the following, replacing ''xx'' with the version of boost downloaded:
Note: the following will overwrite any version of boost currently installed.  It is recommended you either 
remove any version of boost you may have (none on a fresh build) or change --prefix= lines to a different 
directory

      $ tar xjf boost_xx.bz2
      $ cd boost_xx
      $ ./bootstrap.sh --prefix=/usr
      $ sudo ./b2 --without-python --without-mpi --prefix=/usr install

  For the gui:

      $ sudo apt-get install libbz2-dev libqwt5-qt4-dev libqt4-opengl-dev qt4-dev-tools
-# Create the file /usr/share/qt4/mkspecs/features/qwt.prf with the following contents:

        QwtVersion = 5.2.1
        QwtBase = /usr
        LIBS += -L$${QwtBase}/lib -lqwt-qt4
        INCLUDEPATH += $${QwtBase}/include/qwt-qt4

Fedora 16
---------

**Note:** Fedora uses qmake-qt4 instead of qmake (symlink could be used)
-# Run the following

       $ sudo yum install boost-devel gcc-c++ subversion

  For the gui:

      $ sudo yum install sudo yum install qt-devel qwt-devel
-# Create the file /usr/lib/qt4/mkspecs/features/qwt.prf with the following contents:

        QwtVersion = 5.2.2
        QwtBase = /usr

        LIBS += -L$${QwtBase}/lib -lqwt
        INCLUDEPATH += $${QwtBase}/include/qwt


Our code
--------

    $ svn checkout svn://innoculant.micro.med.umich.edu/dev/gr2d/GR-ABM-ODE path/to/my/repo
    $ cd path/to/my/repo
    $ cd simulation
    $ make
    $ ./gr --help

For the gui version, type the following:

    $ cd path/to/my/repo
    $ qmake
    $ make
    $ ./grviz-lung --help
