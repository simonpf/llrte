Tensors, arrays and data
========================

The :code:`llrte/data.h` header provides two template classes to store data:
The :class:`Tensor` class for generic rank-k tensors and the :class:`Array`
class for 1-dimensional data.

Tensors
-------

The Tensor class template implements rank-k tensors.

.. doxygenclass:: llrte::Tensor
   :members:

Data
----

.. doxygenclass:: llrte::Data
   :members:
