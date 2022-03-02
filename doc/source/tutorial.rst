How to use CsstMscData
======================

*CsstMscData* is inherited from *CsstData*.

Example:

.. code-block:: python
    :linenos:

    # set file path
    fp = "MSC_CLB_210525120000_100000000_06_raw.fits"
    # import CsstMscImgData
    from csst.msc import CsstMscImgData
    # read data
    data = CsstMscImgData.read(fp)

To show the info of data:

.. code-block:: python
    :linenos:

    # print info
    print("data: ", data)
    print("instrument: ", data.get_l0keyword("pri", "INSTRUME"))
    print("object: ", data.get_l0keyword("pri", "OBJECT"))

The output:

.. code-block::
    :linenos:

    data:  <CsstMscImgData: MSC CCD16>
    instrument:  MSC
    object:  100000279

Full code is below


.. literalinclude:: ../../examples/example_csstmscdata.py

