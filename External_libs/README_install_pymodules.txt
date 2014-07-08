#############################################
# INSTALLATION of Python modules locally
# (for instance on lxplus at CERN)
#############################################


#########
# NUMPY
#########

tar -xzvf numpy-1.7.0.tar.gz
cd numpy-1.7.0

# if necessary clean out
rm -rf build

# then (to install it in External_libs/numpy-install)

python setup.py build --fcompiler=gnu95
python setup.py install --prefix ../numpy-install

#########
# SCIPY
#########

tar -xzvf scipy-0.12.0.tar.gz
cd scipy-0.12.0

# change setup.py: add 2 lines (at line 145, in function setup_package, before 
# "from numpy.distutils.core import setup")
# (this is to point to the numpy installation directory you just created)
import sys
sys.path.insert(1,'../numpy-install/lib64/python2.6/site-packages');
# NOTE: in the path above you might need to change "python2.6" with your 
# actual python version

# if necessary clean out
rm -rf build

# then (to install it in External_libs/scipy-install)
python setup.py build --fcompiler=gnu95
python setup.py install --prefix=../scipy-install


########################################################
# MATPLOTLIB without GUI (so without interactive plots)
########################################################

# this is a minimal matplotlib, without GUI.
# USE THIS ONLY IF YOU DON'T CARE ABOUT INTERACTIVE PLOTS
# AND ONLY WANT TO OUTPUT PLOTS IN FILES (.png and .eps).
# Otherwise skip this and go to the next.

tar -xzvf matplotlib-1.2.1.tar.gz
cd matplotlib-1.2.1

# change setupext.py: add 1 line (at line 67, after "import sys") (to point to correct numpy installation directory)
sys.path.insert(1,'../numpy-install/lib64/python2.6/site-packages');
# NOTE: in the path above you might need to change "python2.6" with your 
# actual python version

# if necessary clean out
rm -rf build

# create setup.cfg file (to avoid using Tkinter)
cp setup.cfg.template setup.cfg
# modify setup.cfg: uncomment the following lines (to avoid using GTK or Tkinter)
gtk = False
gtkagg = False
tkagg = False
macosx = False


# then (to install it in External_libs/matplotlib-install)
python setup.py build
python setup.py install --prefix=../matplotlib-install


########################################################
# MATPLOTLIB with GUI (so with interactive plots)
########################################################

# USE THIS IF YOU WANT INTERACTIVE PLOTS

# First, install GTK, PyCairo, PyGObject and PyGTK (http://www.pygtk.org/)

# Then:

tar -xzvf matplotlib-1.2.1.tar.gz
cd matplotlib-1.2.1

# change setupext.py: add 1 line (at line 67, after "import sys") (to point to correct numpy installation directory)
sys.path.insert(1,'../numpy-install/lib64/python2.6/site-packages');
# NOTE: in the path above you might need to change "python2.6" with your 
# actual python version

# if necessary clean out
rm -rf build

# then (to install it in External_libs/matplotlib-install)
python setup.py build
python setup.py install --prefix=../matplotlib-install


#################################
# Final .bashrc configuration
#################################

# NOTE: WORKS ONLY IF YOU USE bash SHELL !


# launch the script script_configure_bashrc.sh
./script_configure_bashrc.sh
# and change the paths provided by default if needed (path to 'site-packages'
# folders of numpy, scipy and matplotlib installation directories you just created)

# this will add to your ~/.bashrc some environment variables needed by the python programs 
# and scripts in ../PYTHON_codes_and_scripts/

# finally do
source ~/.bashrc
# (or logout and log in again)
