
OSX Build Instructions  {#readme-osx}
======================

Xcode
-----

-# Update OSX (reboot required)
  -# Click on the apple icon in the upper left corner, then click "Software Update"
  -# Follow the instructions and reboot.
-# Go to http://connect.apple.com/
  -# Sign in with your apple ID (or register for one if you don't)
  -# Download Xcode for your version (Lion requires at least Xcode 4.3.2)
  -# Download Xcode Command-Line Utilities (not required for Snow Leopard or older)
-# Double-click downloaded xcode dmg file
-# Drag and drop xcode.app to Applications folder
-# Double-click downloaded xcode command-line utilities dmg file
-# Double-click dpkg file and follow installation instructions (fairly straight-forward)

Macports
--------

-# Go to http://www.macports.org/install.php
  -# Download dpkg file for your version
-# Double-click and follow installation instructions
-# Open Terminal
  -# Type the following:

            $ sudo xcode-select -switch /Applications/Xcode.app/Contents/Developer
            $ sudo port -v selfupdate
-# To install required libraries for the non-gui version, try these commands:

        $ sudo port install -s boost
-# For the gui version, qt4 and qwt is required:

        $ sudo port install -s qt4-mac qwt-devel
-# All libraries should be installed to /opt/local/lib and header files should be in /opt/local/include

Our code
--------

    $ svn checkout svn://innoculant.micro.med.umich.edu/dev/gr2d/GR-ABM-ODE path/to/my/repo
    $ cd path/to/my/repo
    $ cd simulation
    $ export BOOST_PREFIX=/opt/local
    $ make
    $ ./gr --help

For the gui version, type the following:

    $ cd path/to/my/repo
    $ export BOOST_PREFIX=/opt/local
    $ qmake
    $ make
    $ ./grviz-lung --help

VirtualBox iAtkos Lion
----------------------

**Skip this section if you're not setting up a virtual box for mac osx**

-# Go to http://iatkos.me/ and download iatkos L2 (NOT L2M) dmg
-# Create a new machine with a disk image, set for MacOSX (either 64-bit or 32, didn't seem to make a difference)
-# Open settings of machine
  -# Uncheck "Enable EFI" **VERY IMPORTANT**
  -# Increase base memory (more the better, but remember your host OS needs plenty as well)
  -# Increase number of processors ( > 2 is best)
  -# Increase video ram (at least 64MB, more the better)
  -# In the storage tab, mount the iAtkos dmg file as the CD
-# Run the machine
-#: If you get an error about VT-x / AMD-v, it means you don't have virtualization turned on in the bios and can't vm more than 1 core.  Look in your bios manual for more info
-# Run the standard mac installation process.
   + *Restarting does not work automatically and will crash the vm.  Just reset the machine manually when this happens.*
   + *Do not add your apple ID - this causes a very long sync period on first start up.  You can always do it later and this seems to work better.*
-# Open System Preferences -> Dock
  -# Disable "animate opening applications" <- **Important: gui issues otherwise**
  -# Change "Minimize windows using" to Scale effect <- **Performance reasons**

