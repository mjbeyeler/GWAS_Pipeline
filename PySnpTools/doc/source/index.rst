################################
:mod:`pysnptools` Documentation
################################

PySnpTools: A library for reading and manipulating genetic data.

:synopsis:

* :mod:`.snpreader`: Efficiently read genetic PLINK formats including \*.bed/bim/fam and phenotype files. Also, efficiently read *parts* of files and standardize data.

* :mod:`.kernelreader`: Efficiently create, read, and manipulate kernel data.

* :mod:`.util`: In one line, intersect and re-order IIDs from :mod:`.snpreader`, :mod:`.kernelreader` and other sources. Also, efficiently extract a submatrix from an ndarray.

* :class:`.util.IntRangeSet`: Efficiently manipulate ranges of integers -- for example, genetic position -- with set operators including
  union, intersection, and set difference. 

* :mod:`.pstreader`: Generalizes :mod:`.snpreader` and :mod:`.kernelreader` (provides the efficiency of numpy arrays with some of the flexibility of pandas)

* :mod:`.standardizer`: Specify standardizers for :mod:`.snpreader`.

* :mod:`.kernelstandardizer`: Specify standardizers for :mod:`.kernelreader`.

:Tutorial:

*From PyData Conference, Seattle, 2015*

* `Slides <https://onedrive.live.com/view.aspx?cid=b89ee402873f0f4a&page=view&resid=B89EE402873F0F4A!209845&parId=B89EE402873F0F4A!209839&authkey=!AMAMALlp0TiNuDg&app=PowerPoint&wacqt=undefined>`_
* `Video <https://www.youtube.com/watch?v=KPI6479ctAQ>`_
* `Notebook <https://github.com/MicrosoftGenomics/PySnpTools/blob/master/doc/ipynb/tutorial.ipynb>`_
* `Data file <http://research.microsoft.com/en-us/um/redmond/projects/MSCompBio/PySnpTools/psttutorial.zip>`_

.. automodule:: pysnptools
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:

***********************
:mod:`snpreader` Module
***********************

.. automodule:: pysnptools.snpreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`snpreader.SnpReader`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`snpreader.Bed`
++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Bed
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.SnpData`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpData
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.Pheno`
++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Pheno
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.Ped`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Ped
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.Dat`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Dat
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.Dense`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.Dense
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row

:class:`snpreader.SnpHdf5`
+++++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpHdf5
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`snpreader.SnpNpz`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.snpreader.SnpNpz
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

****************************
:mod:`kernelreader` Module
****************************
.. automodule:: pysnptools.kernelreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`kernelreader.KernelReader`
+++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs, col_property, row_property


:class:`kernelreader.KernelData`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelData
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`kernelreader.KernelNpz`
++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelNpz
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`kernelreader.KernelHdf5`
++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.KernelHdf5
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`kernelreader.Identity`
+++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.Identity
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property

:class:`kernelreader.SnpKernel`
+++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelreader.SnpKernel
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


***********************
:mod:`util` Module
***********************

:mod:`util`
++++++++++++++++++++++++++
.. automodule:: pysnptools.util
    :members:
    :undoc-members:
	:show-inheritance:

:class:`util.IntRangeSet`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.util.IntRangeSet
  :members:
  :undoc-members:
  :show-inheritance:
  :special-members:
  :exclude-members: __and__, __weakref__,__module__,__dict__, __add__

:mod:`util.pheno`
++++++++++++++++++++++++++
.. automodule:: pysnptools.util.pheno
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:

***************************
:mod:`standardizer` Module
***************************

.. automodule:: pysnptools.standardizer
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`standardizer.Standardizer`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Standardizer
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

:class:`standardizer.Unit`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Unit
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

:class:`standardizer.Beta`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Beta
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

:class:`standardizer.Identity`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.Identity
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize,is_constant

:class:`standardizer.DiagKtoN`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.DiagKtoN
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

:class:`standardizer.UnitTrained`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.UnitTrained
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize,is_constant

:class:`standardizer.DiagKtoNTrained`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.standardizer.DiagKtoNTrained
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize,is_constant


******************************************************
:mod:`kernelstandardizer` Module
******************************************************

:class:`kernelstandardizer.KernelStandardizer`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelstandardizer.KernelStandardizer
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs

:class:`kernelstandardizer.DiagKtoN`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelstandardizer.DiagKtoN
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize

:class:`kernelstandardizer.DiagKtoNTrained`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelstandardizer.DiagKtoNTrained
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize,is_constant

:class:`kernelstandardizer.Identity`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.kernelstandardizer.Identity
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs,standardize


***********************
:mod:`pstreader` Module
***********************

.. automodule:: pysnptools.pstreader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:


:class:`pstreader.PstReader`
+++++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstReader
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members: __getitem__
    :exclude-members: copyinputs


:class:`pstreader.PstData`
++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstData
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`pstreader.PstHdf5`
+++++++++++++++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstHdf5
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


:class:`pstreader.PstNpz`
+++++++++++++++++++++++++
.. autoclass:: pysnptools.pstreader.PstNpz
    :members:
    :undoc-members:
	:show-inheritance:
	:special-members:
    :exclude-members: copyinputs, col, col_property, row, row_property


.. only:: html 

***********************
Indices and Tables
***********************

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
