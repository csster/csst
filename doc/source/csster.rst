Guidelines to developers
========================

docstring
---------

We adopt the *numpy* style docstring (see: https://numpydoc.readthedocs.io/en/latest/format.html)

An example: the docstring of numpy.cos()

.. code-block:: python

    def cos(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature, extobj]):
        """
        Cosine element-wise.

        Parameters
        ----------
        x : array_like
            Input array in radians.
        out : ndarray, None, or tuple of ndarray and None, optional
            A location into which the result is stored. If provided, it must have
            a shape that the inputs broadcast to. If not provided or None,
            a freshly-allocated array is returned. A tuple (possible only as a
            keyword argument) must have length equal to the number of outputs.
        where : array_like, optional
            This condition is broadcast over the input. At locations where the
            condition is True, the `out` array will be set to the ufunc result.
            Elsewhere, the `out` array will retain its original value.
            Note that if an uninitialized `out` array is created via the default
            ``out=None``, locations within it where the condition is False will
            remain uninitialized.
        **kwargs
            For other keyword-only arguments, see the
            :ref:`ufunc docs <ufuncs.kwargs>`.

        Returns
        -------
        y : ndarray
            The corresponding cosine values.
            This is a scalar if `x` is a scalar.

        Notes
        -----
        If `out` is provided, the function writes the result into it,
        and returns a reference to `out`.  (See Examples)
        References
        ----------
        M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions.
        New York, NY: Dover, 1972.

        Examples
        --------
         np.cos(np.array([0, np.pi/2, np.pi]))
        array([  1.00000000e+00,   6.12303177e-17,  -1.00000000e+00])
         # Example of providing the optional output parameter
         out1 = np.array([0], dtype='d')
         out2 = np.cos([0.1], out1)
         out2 is out1
        True
         # Example of ValueError due to provision of shape mis-matched `out`
         np.cos(np.zeros((3,3)),np.zeros((2,2)))
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        ValueError: operands could not be broadcast together with shapes (3,3) (2,2)
        Class docstring:
        Functions that operate element by element on whole arrays.
        To see the documentation for a specific ufunc, use `info`.  For
        example, ``np.info(np.sin)``.  Because ufuncs are written in C
        (for speed) and linked into Python with NumPy's ufunc facility,
        Python's help() function finds this page whenever help() is called
        on a ufunc.
        A detailed explanation of ufuncs can be found in the docs for :ref:`ufuncs`.
        **Calling ufuncs:** ``op(*x[, out], where=True, **kwargs)``
        Apply `op` to the arguments `*x` elementwise, broadcasting the arguments.
        The broadcasting rules are:
        * Dimensions of length 1 may be prepended to either array.
        * Arrays may be repeated along dimensions of length 1.

        Parameters
        ----------
        *x : array_like
            Input arrays.
        out : ndarray, None, or tuple of ndarray and None, optional
            Alternate array object(s) in which to put the result; if provided, it
            must have a shape that the inputs broadcast to. A tuple of arrays
            (possible only as a keyword argument) must have length equal to the
            number of outputs; use None for uninitialized outputs to be
            allocated by the ufunc.
        where : array_like, optional
            This condition is broadcast over the input. At locations where the
            condition is True, the `out` array will be set to the ufunc result.
            Elsewhere, the `out` array will retain its original value.
            Note that if an uninitialized `out` array is created via the default
            ``out=None``, locations within it where the condition is False will
            remain uninitialized.
        **kwargs
            For other keyword-only arguments, see the :ref:`ufunc docs <ufuncs.kwargs>`.

        Returns
        -------
        r : ndarray or tuple of ndarray
            `r` will have the shape that the arrays in `x` broadcast to; if `out` is
            provided, it will be returned. If not, `r` will be allocated and
            may contain uninitialized values. If the function has more than one
            output, then the result will be a tuple of arrays.

        """


TODOs
-----

Revise your code/docstring if you get a "TODO" mark in your code, plz.

`logging`
---------

Use `logging` module to print verbose info.

Environment variables
---------------------

Define some environment variables for testing.
