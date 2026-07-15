Installation
===============================

===============================
Via PIP
===============================

| You need a working git and python3 installations for this guide to work.

You can install unduwave simply via pip by typing:

.. code-block:: console

   pip install unduwave

After everything is finished, you can import unduwave via:

.. code-block:: console

   import unduwave as uw

| After this, you can open a terminal in the "scripts/First-Examples" subfolder - or navigate there from inside a terminal using the "cd" command - and run:

.. code-block:: console

   python3 example_wave.py 

The script should run your first WAVE simulation.

===============================
Via github
===============================

Clone the repository

.. code-block:: console

   git clone https://github.com/SteJSch/unduwave

Navigate to the package folder and build the dist via:

.. code-block:: console

   python -m build
   python -m pip install . 

This should build and install unduwave as a python package on your system.
   

