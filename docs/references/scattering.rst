Scattering
==========

*Improving workflow performance with embarrassingly parallel tasks*

*This page is under construction, please check for more details later*

Basic information:

- Janis only officially supports scattering on one field.

- Scattering implicitly occurs when an array of values is passed to an input that
only accepts single values (of the same type). A merge always occurs in the proceeding
step, though it may be *rescattered*.


WDL
*******
Calls in WDL are explicitly scattered with essentially a for-loop. This calls for
intermediary aliasing of the variable. It attempts to intelligently do this based
on the field name to scatter on where possible.

.. code-block:: wdl

   scatter (b in bams) {
      call G.gatk4sortsam as sortsam {
       input:
         bam=b,
     }
   }


Scattering by more than one field
*********************************

This is currently unsupported.

In CWL this could be achieved natively, in WDL it would probably require creating
nested sub workflows with this logic achieved through the `zip` and `cross` functions,
then you could let the WDL engine do the implicit ordered gathering.
