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

A pipeline should have the structure like below. (see csst/)

.. code-block:: python

    from csst.core.processor import CsstProcessor

    class CsstProcDemo(CsstProcessor):
        def prepare(self, **args):
            # prepare the environment
            # for example, if you make use of some third-party software like SEXTRACTOR,
            # do your preparation here.
            pass

        def run(self, CsstData, *args, **kwargs):
            # run your pipeline here
            # make sure that your input data should be a child class instance of CsstData.
            pass

        def cleanup(self, **kwargs):
            # clean up environment
            pass
