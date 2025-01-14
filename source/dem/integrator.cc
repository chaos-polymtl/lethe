#include <dem/integrator.h>

template class Integrator<2, DEM::DEMProperties::PropertiesIndex>;
template class Integrator<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class Integrator<3, DEM::DEMProperties::PropertiesIndex>;
template class Integrator<3, DEM::CFDDEMProperties::PropertiesIndex>;
