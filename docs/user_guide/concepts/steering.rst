.. _concept_steering:

Steering: Beam Geometry
========================

The *steering information* (``stf``) describes the physical geometry of all treatment beams:
their directions, the positions of individual pencil beams / spots, and — for particle therapy —
the energy layers within each beam.

Generating the steering information
------------------------------------

:func:`~pyRadPlan.generate_stf` reads beam directions and other settings from the plan's
``prop_stf`` dictionary and produces a fully validated
:class:`~pyRadPlan.stf.SteeringInformation` object:

.. code-block:: python

    from pyRadPlan import generate_stf, IonPlan

    pln = IonPlan(radiation_mode="protons", machine="Generic")
    pln.prop_stf = {
        "gantry_angles": [0, 90, 270],   # degrees; one beam per angle
        "couch_angles":  [0,  0,   0],
        "bixel_width":   5,              # mm lateral spot spacing
        "generator":     "IMPT",         # default for ions
    }

    stf = generate_stf(ct, cst, pln)

Available generators
--------------------

The generator is selected via ``prop_stf["generator"]`` using its short name. If omitted, a
default suitable for the plan's radiation mode is chosen. The generators compatible with a
given plan can be listed with :func:`~pyRadPlan.stf.generators.get_available_generators`.

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Short name
     - Radiation modes
     - Description
   * - ``IMPT``
     - protons, helium, carbon, oxygen
     - :class:`~pyRadPlan.stf.generators.StfGeneratorIMPT` — Intensity-modulated particle
       therapy. Places a lateral spot grid (``bixel_width`` / ``lateral_spot_spacing``) over the
       target projection and selects, per ray, all energy layers whose Bragg peaks fall inside
       the target. Energy layers are thinned to ``longitudinal_spot_spacing`` (default 2 mm).
   * - ``ionSingleSpot``
     - protons, helium, carbon, oxygen
     - :class:`~pyRadPlan.stf.generators.StfGeneratorIonSingleSpot` — A single ray through the
       isocenter with one spot at a fixed energy (``fix_energy``; the closest machine energy is
       used). Useful for pencil-beam studies, tests, and debugging.
   * - ``photonIMRT``
     - photons
     - :class:`~pyRadPlan.stf.generators.StfGeneratorPhotonIMRT` — Ray-bixel geometry for
       photon IMRT. Covers the target projection with a grid of bixels of size ``bixel_width``,
       each at a single nominal energy (default 6 MV).
   * - ``photonOpenFields``
     - photons
     - :class:`~pyRadPlan.stf.generators.StfGeneratorPhotonCollimatedSquareFields` — One central
       ray per beam representing a collimated square open field of width ``field_width``.
   * - ``photonSingleBixel``
     - photons
     - :class:`~pyRadPlan.stf.generators.StfGeneratorPhotonSingleBixel` — A single bixel through
       the isocenter per beam. Optionally supports field-based dose calculation with a custom
       field shape and beam-limiting devices.
   * - ``VHEE``
     - VHEE
     - :class:`~pyRadPlan.stf.generators.StfGeneratorVHEE` — Very-high-energy electrons. Uses a
       ray-bixel grid over the target like IMRT, but with a single fixed energy per beam (default
       200 MeV).

Custom generators can be added by subclassing
:class:`~pyRadPlan.stf.generators.StfGeneratorBase` and registering them with
:func:`~pyRadPlan.stf.generators.register_generator`.

Object hierarchy
-----------------

The steering information is organized as a three-level hierarchy:

.. code-block:: text

    SteeringInformation
    └── Beam  (one per gantry/couch angle pair)
        └── Ray  (one per lateral pencil-beam position)
            └── Beamlet  (one per energy / focus setting)

SteeringInformation
~~~~~~~~~~~~~~~~~~~

:class:`~pyRadPlan.stf.SteeringInformation` is the top-level container.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``beams``
     - List of :class:`~pyRadPlan.stf.Beam` objects.
   * - ``num_of_beams``
     - Total number of beams.
   * - ``num_of_rays``
     - Total ray count across all beams.
   * - ``total_number_of_bixels``
     - Total beamlet count. This is the column dimension of the dose influence matrix.
   * - ``bixel_beam_index_map``
     - Integer array mapping every global bixel index to its beam index.

Beam
~~~~

:class:`~pyRadPlan.stf.Beam` represents a single irradiation direction.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``gantry_angle``
     - Gantry rotation in degrees (IEC convention).
   * - ``couch_angle``
     - Couch rotation in degrees.
   * - ``bixel_width``
     - Lateral spot spacing in mm.
   * - ``sad``
     - Source-to-axis distance in mm.
   * - ``iso_center``
     - Isocenter coordinates ``[x, y, z]`` in world coordinates.
   * - ``source_point``
     - Source position in world coordinates.
   * - ``rays``
     - List of :class:`~pyRadPlan.stf.Ray` objects.

Ray
~~~

:class:`~pyRadPlan.stf.Ray` represents one lateral pencil-beam direction within a beam.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``ray_pos``
     - Ray position in world coordinates (mm).
   * - ``ray_pos_bev``
     - Ray position in beam's-eye-view (BEV) coordinates.
   * - ``beamlets``
     - List of :class:`~pyRadPlan.stf.Beamlet` objects (one per energy layer).

Beamlet
~~~~~~~

:class:`~pyRadPlan.stf.Beamlet` is the smallest unit — a single pencil beam at a fixed
energy.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``energy``
     - Beam energy in MeV.
   * - ``weight``
     - Initial fluence weight (updated by the optimizer).
   * - ``num_particles_per_mu``
     - Conversion factor from monitor units to particle count.
   * - ``min_mu`` / ``max_mu``
     - Monitor-unit bounds for the optimizer.
   * - ``range_shifter``
     - Optional range-shifting device description.

Coordinate systems
------------------

pyRadPlan uses the **LPS** (Left-Posterior-Superior) coordinate system for world coordinates,
consistent with DICOM and SimpleITK. Beam's-eye-view (BEV) coordinates follow the IEC 61217
convention.

matRad interoperability
-----------------------

.. code-block:: python

    matrad_stf = stf.to_matrad()  # list of dicts, one per beam
