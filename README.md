# PGSuper IFC Extensions
Experimental features for PGSuper/PGSplice

This project is a sandbox for developing experimental IFC features for [BridgeLink](github.com/wsdot/bridgelink) [PGSuper and PGSplice](github.com/wsdot/pgsuper). 

## Things to remember for the IDS

* Beams must be typed with IfcBeamType
* Beams must have Pset_PrecastConcreteELementGeneral property set with
    * ReleaseStrength
    * (final strength?)
    * DesignLocationNumber and its expected to be "Span i, Girder j" so use a reg expression
* Beams must be classified with usBridge_GirderPrestressedConcrete from US Bridge data dictionary on bSDD
* Beams must have Pset_ConcreteElementGeneral with
    * AssemblyPlace = FACTORY
    * CastingMethod = PRECAST


We want to have an input IDS to validate models before importing them.

We also want to output an IDS that contains things like the design values for concrete strength, strand forces, camber, etc so the proper design results can be evaluated against the final model which will be created in a model authoring tool

