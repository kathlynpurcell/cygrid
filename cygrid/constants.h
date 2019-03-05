#ifndef CONSTANTS_H
#define CONSTANTS_H

// --------------------------------------------------------------------
//
// title                  :constants.h
// description            :Some constants to be used in Cygrid.
// author                 :Benjamin Winkel, Lars Flöer & Daniel Lenz
//
// --------------------------------------------------------------------
//  Copyright (C) 2010+ by Benjamin Winkel, Lars Flöer & Daniel Lenz
//  bwinkel@mpifr.de, mail@lfloeer.de, dlenz.bonn@gmail.com
//  This file is part of cygrid.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// --------------------------------------------------------------------

// Define maximal y-size of pixel indices
// this is necessary to have a quick'n'dirty hash for the xpix-ypix pairs
// otherwise we would need to provide a hash function to unordered_map,
// (and we don't want to use the true y-size of the map, to avoid a dynamical
// variable for this)
#define MAX_Y (1073741824)  // 2^30

// Constants needed for HPX
#define NESTED 1
#define RING 2

#define PI (3.1415926535897932384626433832)
#define TWOTHIRD (2.0 / 3.0)
#define TWOPI (2 * PI)
#define HALFPI (PI / 2.)
#define INV_TWOPI (1.0 / TWOPI)
#define INV_HALFPI (1. / HALFPI)
#define DEG2RAD (PI / 180.)
#define RAD2DEG (180. / PI)

#endif

