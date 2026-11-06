.. _concept_dose_calculation:

Dose Calculation
================

Dose calculation in pyRadPlan produces a *dose influence matrix* (``dij``) that encodes the
dose deposited in every voxel by every beamlet. The influence matrix is the central data
structure for optimization: the dose for any fluence vector ``w`` is simply the matrix–vector
product ``dij · w``.

Two top-level functions are provided:

- :func:`~pyRadPlan.calc_dose_influence` — compute the full dose influence matrix, storing the
  dose of every beamlet as a separate column.
- :func:`~pyRadPlan.calc_dose_forward` — compute the dose of one given fluence vector directly,
  accumulating the weighted beamlet doses instead of storing them.

The forward calculation saves memory, since no influence matrix is stored. It is useful for recalculating a plan, e.g.
with a different engine or dose grid, or when the influence matrix would not fit into memory.
To evaluate the dose of an optimized fluence with an existing ``dij``, use
``dij.compute_result_ct_grid(fluence)`` instead (see below).

.. code-block:: python

    from pyRadPlan import calc_dose_influence, calc_dose_forward

    # Full influence matrix (used for optimization)
    dij = calc_dose_influence(ct, cst, stf, pln)

    # Direct dose of a given fluence; returns result images on the CT grid
    result = calc_dose_forward(ct, cst, stf, pln, weights=fluence)
    physical_dose_ct = result["physical_dose"]

Both functions select the appropriate engine automatically from ``pln.prop_dose_calc["engine"]``.

The Dij object
--------------

:class:`~pyRadPlan.dij.Dij` stores the dose influence matrix together with the grid metadata
needed to relate dose voxels back to CT voxels.

Dose matrices
~~~~~~~~~~~~~

The primary attribute is ``physical_dose``. Influence matrices are stored in object-array
containers so that single- and multi-scenario calculations use the same layout; each contained
matrix has shape ``(num_voxels, num_bixels)`` and is commonly stored in CSC sparse format.
Additional matrices are available for biological dose calculation in particle therapy:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Attribute
     - Description
   * - ``physical_dose``
     - Physical dose contribution per voxel per beamlet (Gy / unit fluence).
   * - ``let_dose``
     - LET-weighted dose matrix for LET-based optimization (keV/µm · Gy / unit fluence).
   * - ``alpha_dose``
     - Per-beamlet α-dose for the linear-quadratic (LQ) biological model (α in Gy⁻¹ times
       dose, i.e. 1 / unit fluence).
   * - ``sqrt_beta_dose``
     - Per-beamlet √β-dose for the LQ model (√β in Gy⁻¹ times dose, i.e. 1 / unit fluence).

The set of available matrices determines which :ref:`quantities <concept_quantities>` can be
resolved during optimization. Every engine, photon engines included, stores the biological model it used as
``dij.bio_model``; for a constant-RBE model no matrices are needed and the constant is
available as ``dij.rbe``. An engine that cannot compute a model's declared outputs (for
example a photon engine with an LQ model) still records the model but logs a warning that
only physical dose is available. The model decides how ``rbe_x_dose`` is derived (LQ inversion of
the effect, or ``rbe * physical_dose``).

Biological models and LET
-------------------------

Which of the additional matrices a particle engine computes follows from ``pln.bio_model``
(see :ref:`concept_bio_model`) and the machine data:

* Supported biological influence matrices are computed when the evaluator declares them. The
  current ``Dij`` schema supports the ``alpha_dose`` / ``sqrt_beta_dose`` declared by the in-tree
  LQ evaluators.
* ``let_dose`` is computed when the machine provides LET kernels (pencil beam) or the biological
  model requires LET. How the LET is obtained depends on the engine: pencil-beam engines
  interpolate it from the machine's LET kernels and therefore need them, while Monte Carlo
  engines such as FRED can score the LET themselves during the simulation. An LET-based model can
  therefore be used with FRED even if the machine carries no LET kernels.

Both can be overridden through the engine parameters ``calc_bio_dose`` and ``calc_let``,
which accept ``"auto"`` (default), ``True`` or ``False``:

.. code-block:: python

    pln.bio_model = "MCN"
    # keep only physical dose and LET; evaluate the LET-based model later
    pln.prop_dose_calc = {"engine": "HongPB", "calc_bio_dose": False}

Setting ``calc_bio_dose=True`` without an evaluator providing influence quantities raises an
error. An evaluator declaring a quantity that ``Dij`` cannot currently store also fails during
setup. Switching ``calc_let`` off does not affect LET-based models, which still receive the LET
kernels internally. The FRED engine derives biological matrices from its scored LET and therefore
supports LET-based models only.

During the calculation the model is bound to the machine and the tissue properties of the
voxels (from the structure set, with ``dij.alphax`` / ``dij.betax`` holding their reference
parameters) through a :class:`~pyRadPlan.biology.models.BioModelEvaluator`. Kernel- and
table-based models carry base data for a finite set of tissue classes, and select the class
of each voxel by exact match against them
(:class:`~pyRadPlan.biology.models.ExactClassLookup`). What is compared are the fields
identifying a tissue class of the model, not always the reference parameters: LEM base data
is tabulated per ``(alpha_x, beta_x, d_t, r_nucleus, version)``, while the MKM z\* kernels
only depend on the ``(r_nucleus, r_domain)`` geometry — its ``alpha_0`` / ``beta_0`` are taken
from the tissue per voxel. Only the tissue classes of the kernel set a model resolves to are
considered, and a tissue of another provenance (e.g. MKM parameters for LEM I kernels) is
rejected, so a machine carrying classes for several models cannot serve the base data of
one through the tissue class of another. The engine runs
:meth:`~pyRadPlan.biology.models.BiologicalModel.check_tissues` before the calculation. LEM
kernels identify their rows by the reference ``alpha_x`` / ``beta_x`` stored with them; the
machine's declared tissue classes are matched onto these rows by value, and a machine without
declared classes uses the kernels' own. See :ref:`voi_tissue_properties` for setting the
structures up and :ref:`biological_base_data` for the kernel sets and their metadata.

A model that cannot be applied to a calculation — the machine lacks its kernel quantities, a
voxel's tissue matches no class of the base data, or the reference parameters are outside the
model's domain — only fails the calculation when biological dose is requested explicitly
(``calc_bio_dose=True``). With ``"auto"`` or ``False`` a warning names the reason and the
physical dose influence is computed alone: the dij then records no model
(``dij.bio_model is None``) and carries no biological matrices, so ``rbe_x_dose`` is reported
as unavailable and the legacy objective conversion keeps physical-dose objectives on physical
dose. This keeps physical-dose workflows working with the per-modality default models.
Invalid reference coefficients (negative, infinite or NaN ``alpha_x`` / ``beta_x``) always
raise.

Grid information
~~~~~~~~~~~~~~~~

Dose is typically computed on a coarser *dose grid* (for speed) and then interpolated back to
the *CT grid* for display:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Attribute
     - Description
   * - ``dose_grid``
     - :class:`~pyRadPlan.core.Grid` describing the voxel grid on which dose was computed.
   * - ``ct_grid``
     - :class:`~pyRadPlan.core.Grid` of the original CT image.
   * - ``num_of_voxels``
     - Number of voxels in the dose grid (= ``dose_grid.num_voxels``).
   * - ``total_num_of_bixels``
     - Number of beamlets (= columns of the dose matrix).

Reconstructing dose on the CT grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Compute dose for optimized fluence and resample to CT grid
    result = dij.compute_result_ct_grid(fluence)
    physical_dose_ct = result["physical_dose"]  # SimpleITK image on the CT grid

    # Inspect which quantities are available
    print(dij.quantities)   # e.g. ['physical_dose', 'alpha_dose', 'sqrt_beta_dose']

Dose engines
------------

The dose engine is selected via ``pln.prop_dose_calc["engine"]``. All engines inherit from a
common abstract base and implement the same interface, so switching between them requires only
one line:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Engine key
     - Description
   * - ``"SVDPB"``
     - SVD pencil beam for photons.
   * - ``"HongPB"``
     - Hong's analytical pencil beam for protons and ions.
   * - ``"FRED"``, ``"TOPAS"``
     - Monte Carlo engines interfacing external codes (installation required):

       - ``"FRED"`` — `FRED <https://www.fred-mc.org>`_. Validated for protons only; ion
         calculations warn that their doses do not yet reproduce the reference.
       - ``"TOPAS"`` — `TOPAS <https://opentopas.readthedocs.io>`_ for photons, VHEE, protons and
         ions.

.. code-block:: python

    pln.prop_dose_calc = {"engine": "HongPB"}
    dij = calc_dose_influence(ct, cst, stf, pln)

Pencil-beam engines are fast and deterministic and are the usual choice for computing the
influence matrix for optimization. Monte Carlo engines model particle transport in
heterogeneous tissue more accurately at a higher computational cost. All engine parameters are listed in the API reference of
the respective engine class.

SVD pencil beam
~~~~~~~~~~~~~~~

:class:`~pyRadPlan.dose.engines.PhotonPencilBeamSVDEngine` follows the singular value
decomposition method of Bortfeld, Schlegel & Rhein (1993). The depth-dependent photon
pencil-beam kernel is decomposed into a few depth-independent lateral kernels with
depth-dependent weights. The radiological depth of each voxel is obtained by ray tracing
through the CT, the bixel fluence is convolved with the lateral kernels in the beam's eye
view, and the result is evaluated at the voxel positions and depths. Low-dose entries far
from a beamlet can be thinned by random sampling (*dij sampling*) to keep the matrix sparse.

Hong pencil beam
~~~~~~~~~~~~~~~~

:class:`~pyRadPlan.dose.engines.ParticleHongPencilBeamEngine` implements the particle
pencil-beam algorithm of Hong et al. (1996). Each spot is modelled as a tabulated integral
depth dose, interpolated at the water-equivalent depth of the voxel, multiplied by a lateral
Gaussian whose width grows with depth through multiple Coulomb scattering and is broadened by
the initial spot size. LET and biological kernels of the machine are
interpolated the same way. Heterogeneities are only taken into account along the central ray
of each spot.

Monte Carlo engines
~~~~~~~~~~~~~~~~~~~

The Monte Carlo engines transport individual particles through the voxelized patient and
therefore model lateral scattering in heterogeneous tissue and nuclear interactions
explicitly. pyRadPlan does not simulate the particles itself but interfaces external codes,
which must be installed separately: it writes the CT with a HU-to-material conversion and the
plan's fields, energy layers and spots as input for the code, and reads the scored dose back
into a ``Dij``. Both engines compute a single (nominal) scenario and score the dose-averaged
LET themselves, so LET-based biological models work without LET kernels in the machine.

- **FRED** — :class:`~pyRadPlan.dose.engines.ParticleFredMCEngine` interfaces the
  GPU-accelerated Monte Carlo code FRED (Schiavi et al. 2017), which pyRadPlan runs directly.
  For an influence matrix FRED simulates ``num_histories_per_beamlet`` primaries per spot; a
  forward calculation simulates ``num_histories_direct`` primaries for the whole plan.
- **TOPAS** — :class:`~pyRadPlan.dose.engines.ParticleTOPASMCEngine` interfaces TOPAS, a
  Geant4-based Monte Carlo toolkit (Perl et al. 2012; Faddegon et al. 2020). As the TOPAS monte carlo
  code is not optimized for speed, it is meant for forward calculation of a given plan. Therefore the calculation
  runs in two steps: with ``external_calculation=True`` pyRadPlan writes the input to
  ``simu_dir``, sampling ``num_histories_direct`` primaries from the spot weights and splitting
  them into ``num_runs`` independent runs. After running TOPAS, setting
  ``external_calculation`` to the output folder reads the scored dose, sums the runs per beam
  and stores the run-to-run variance. The result holds one column per beam.
   TOPAS has the option to score proton LET.

Backend-agnostic computation
-----------------------------

Dose engine and quantity code use the
`Python Array API standard <https://data-apis.org/array-api/latest/>`_ where practical. The
preferred compute backend is the ``xp`` sub-configuration of the global
:data:`pyRadPlan.settings` and can be set via environment variables
(``PYRADPLAN_XP_PREFER_GPU``, ``PYRADPLAN_XP_PREFERRED_GPU_ARRAY_BACKEND``), a ``.env`` file, or
at runtime:

.. code-block:: python

    from pyRadPlan import settings

    settings.xp.prefer_gpu = True
    settings.xp.preferred_gpu_array_backend = "cupy"
    dij = calc_dose_influence(ct, cst, stf, pln)

The ``Dij`` object can also convert its influence matrices explicitly with
``dij.to_namespace("cupy")`` or another supported namespace.

.. _jit-compiled-kernels:

Jit-compiled kernels
~~~~~~~~~~~~~~~~~~~~

A handful of hot-path computations (currently in the Siddon ray tracer) are written once
as plain Array API code and additionally wrapped in
:func:`pyRadPlan.core.xp_utils.jittable`. This decorator gives such a kernel jit-compiled
execution paths per backend — a `numba <https://numba.pydata.org>`_ implementation for
NumPy, or the backend's own tracer (``jax.jit``, ``torch.compile``) for backends where the
kernel is written in a jit-traceable, static-shape style — while still falling back to the
plain implementation everywhere else, and whenever a compiled path is unavailable or fails.

Which backends actually take the jit-compiled path is controlled by
``settings.xp.jit_backends``, a comma-separated list (``PYRADPLAN_XP_JIT_BACKENDS``,
default ``"numpy,jax"``):

.. code-block:: python

    from pyRadPlan import settings

    settings.xp.jit_backends = "numpy,jax,torch"   # also let torch.compile try
    settings.xp.jit_backends = ""                  # disable jitting everywhere

The NumPy path requires the ``[numba]`` extra (``pip install "pyRadPlan[numba]"``); without
it, NumPy silently runs the plain implementation. If a backend's jit fails to compile or
raises on its first execution (this can happen with ``torch.compile`` when the platform's
Triton install is missing or incompatible), pyRadPlan emits a ``RuntimeWarning`` once and
continues with the plain implementation for the remainder of the run — jitting is always a
performance option, never a correctness requirement.

Multi-scenario support
-----------------------

When the plan carries an uncertainty scenario model (``pln.scenarios``, see
:ref:`concept_uncertainty`), ``calc_dose_influence`` resolves it into a
:class:`~pyRadPlan.scenarios.ScenarioSet`, checks the set against the engine's
``supported_uncertainties`` and returns a ``Dij`` whose quantity matrices are flat object-array
containers with one influence matrix per computed scenario, addressed by the scenario's flat
index (``dij.scenarios`` records the set). Scenarios sharing CT phase and setup shift share the
ray tracing; the range perturbation is applied per scenario. The quantities evaluate every
scenario of the set and the optimization problem aggregates over them.
