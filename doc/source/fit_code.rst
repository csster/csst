Fit your code to this project
=============================
Take the main survey (MSC) as an example.
A demo structure is shown below.

.. code-block:: shell

    msc
    ├── __init__.py
    ├── data.py
    │   └── CsstMscImgData class
    └── processor.py
        └── CsstMscImgL0Proc class


The data class
--------------

For each instrument, a specific data class should be constructed.
`CsstData` class is the root class.

.. code-block:: shell

    CsstData
    └── CsstMscData
        └── CsstMscImgData

The pipeline
------------

A pipeline should have the structure like below.

.. code-block:: shell

    CsstMscImgL0Proc
    ├── prepare()
    ├── run()
    └── cleanup()
