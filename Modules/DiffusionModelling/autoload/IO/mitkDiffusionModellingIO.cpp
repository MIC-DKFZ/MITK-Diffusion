/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include <usModuleActivator.h>

// MitkDiffusionModellingIO has no code of its own: MitkCore auto-loads it, and loading it
// loads MitkDiffusionModelling, whose activator registers the readers and writers (see
// CMakeLists.txt). The linker keeps MitkDiffusionModelling as a dependency only if a symbol of
// it is used (--as-needed in MITK's linker flags; Windows); the one used is the entry point of
// that activator, which US_EXPORT_MODULE_ACTIVATOR exports.
extern "C" us::ModuleActivator *_us_module_activator_instance_MitkDiffusionModelling();

us::ModuleActivator *(*mitkDiffusionModellingIODependency)() =
  &_us_module_activator_instance_MitkDiffusionModelling;
