Installation
===============================

===============================
Linux
===============================
| You need working git and python3 installations for thus guide to work.

Open a terminal in the folder where you want to have unduwave. Type:

.. code-block:: console

   git clone --recursive https://github.com/SteJSch/unduwave

After everything is finished, you can go to the unduwave folder and run in a terminal:

.. code-block:: console

   pip3 install -r requirements.txt

| This installs all required packages (consider adding a virtual environment).

| After this, you can go the "scripts" subfolder and run:

.. code-block:: console

   python3 example_wave.py 

The script should run your first WAVE simulation.

If you encounter errors, you may have to recompile WAVE. 
| To do that, first go to the folder "\\bin" inside the "External-Software/WAVE" directory and delete all files. Then, in the terminal, do

.. code-block:: console

   export WAVE_INCL="path/to/unduwave/"

Then type: 

.. code-block:: console

   python3 python/make_wave.py
   
The compilation should now run. After that, go back to the scripts folder and re-run the "example_wave.py" file. 
   
===============================
Windows
===============================

Preliminaries 
----------------------------------

| First, read every step before you start doing things!

| Then, start by installing windows if you havent done that :)

Now, we will begin:

| Download and install Github. For example from here: `Git <https://git-scm.com/downloads/win>`_ .

| Download and install Cygwin. For example from here: `Cygwin <https://cygwin.com/install.html>`_ . 

| Keep the Cygwin installer file. If you need it to install additional packages later, you can just run the installer file to do that. During the installation, choose: "Install from Internet" when asked and choose some mirror that suites you. 

| After the Cygwin installation, the package manager will automatically pop up:

.. figure:: pics/cygwin1.png
   :scale: 50 %
   :alt: map to buried treasure

   Cygwin package manager.

| Under "View" choose "Full". Under "Search" search for the following packages:

.. figure:: pics/cygwin2.png
   :scale: 50 %
   :alt: map to buried treasure

   Cygwin packages to install.

Select the drop-down arrow to the right to select a version of a given program and to add it to the installation list. After you added all packages you can click "next" and the installation will start.

.. figure:: pics/cygwin3.png
   :scale: 50 %
   :alt: map to buried treasure

   How to choose to install a given package.

| The installation may take a long time. Do some sports meanwhile to get your head free :D

| Once the installation is finished, go to the cygwin installation folder and double click on the "cygwin.bat" - Cygwin will be started the first time and does some setting-things. You can close the window afterwards if you like.

| Depending on how you run python programs, you may have to install python again separately. In this guide I will run the scripts from the windows powershell. Search in the windows-search window on the lower left for powershell, start it and type "python3". If python3 is not installed, you will be redirected to the windows store where you can install it. After installation restart the powershell and type again "python3" - the python console should now open.

Github and Cloning
----------------------------------

In the installation folder of Cygwin, you will find a folder "home/user" with user being your user name.
Navigate into this folder with windows-file browser, right click and select "Open Git Bash here"
In the bash write:

.. code-block:: console

   git clone --recursive https://github.com/SteJSch/unduwave

If you do not have the Open Git Bash option on right-click, search for the git-bash in the windows search function and navigate inside of it to the right folder in the Cygwin main folder (using +cd folder1/folder2/.../)

| If this worked, the folder unduwave should appear and all submodules (wave and undumag) should be initialized (in the External-Software folder in unduwave) This may take some time. Get a coffee (or tea).

Test and Compilation
----------------------------------

| Go to the folder "unduwave/External-Software" and open the file "where_is_cygwin_installation.txt" . In the first line, write the path of your cygwin installation (this is a hard hack, I know. But only temporary :) ). Use the format given, with the quotes and the path terminated by a backslash. 

| In the second line, write the path where you put unduwave inside the "home/user/" folder inside cygwin. If unduwave lies at the "cygwin/home/user/" you can leave it empty
   
This will install the required python packages to run the scripts.

| To test the installation, go to the folder "unduwave/" inside the windows powershell (this is what I used to start the python scripts). Start the powershell by searching for it in the windows search field (lower left corner) and then navigate to the correct folder inside the powershell using "cd folder1/folder2/..." command. Then do:

.. code-block:: console

   pip3 install -r requirements.txt
   
This installs all the required python packages to run the scripts.

| Then navigate to the sub-folder scripts and once you are there, do: 

.. code-block:: console

   python3 example_wave.py 

If you are (very) lucky, the program will run, wave will be called and execute normally. 
Otherwise, you will have to compile WAVE and Undumag. 

| To do that, first go to the folder "\\bin" inside the "unduwave\\External-Software\\WAVE" directory and delete all files.

| Then go to the main folder of your cygwin installation, double click on "cygwin.bat" - this opens the cygwin console

| You start at your cygwin home-folder (i.e. in "cygwin/home/user/") and navigate to the folder "unduwave/External-Software/WAVE". In the cygwin console then type: 

.. code-block:: console

   export WAVE_INCL="\path\to\unduwave\External-Software\WAVE"

Where you can use the relative path within you cygwin user directory. So, if unduwave lies in "C:\\cygwin64\\home\\user\\unduwave" you can simply write :

.. code-block:: console

   export WAVE_INCL="unduwave\External-Software\WAVE"

Then type: 

.. code-block:: console

   python3 python/make_wave.py

The compilation should run smoothly, might take a long time. Get yourself whatever you need, a beer or cacao.  
After it is done, go back to the script folder (in the windows powershell) and retry to run the "example_wave.py" file.

If you are lucky, you now have a working unduwave version on you windows PC. If that installation guide was too much for you, perhaps you should consider changing to Linux :)
