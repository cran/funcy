/* MLPACK 0.2
 *
 * Copyright (c) 2008, 2009 Alexander Gray,
 *                          Garry Boyer,
 *                          Ryan Riegel,
 *                          Nikolaos Vasiloglou,
 *                          Dongryeol Lee,
 *                          Chip Mappus, 
 *                          Nishant Mehta,
 *                          Hua Ouyang,
 *                          Parikshit Ram,
 *                          Long Tran,
 *                          Wee Chin Wong
 *
 * Copyright (c) 2008, 2009 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */
/**
 * @file discrete.cc
 *
 * Helpers for discrete math (implementation).
 */

#include "discrete.h"


double math::Factorial(int d) {
  double v = 1;
  
  DEBUG_ASSERT(d >= 0);
  
  for (int i = 2; i <= d; i++) {
    v *= i;
  }
  
  return v;
}

void math::MakeIdentityPermutation(fl__index_t size, fl__index_t *array) {
  for (fl__index_t i = 0; i < size; i++) {
    array[i] = i;
  }
}

void math::MakeRandomPermutation(fl__index_t size, fl__index_t *array) {
  // Regular permutation algorithm.
  // This is cache inefficient for large sizes; large caches might
  // warrant a more sophisticated blocked algorithm.
  
  if (size == 0) {
    return;
  }
  
  array[0] = 0;

  /*
  for (fl__index_t i = 1; i < size; i++) {
    fl__index_t victim = rand() % i;
    array[i] = array[victim];
    array[victim] = i;
  }
  */
}

void math::MakeInversePermutation(fl__index_t size,
    const fl__index_t *original, fl__index_t *reverse) {
  for (fl__index_t i = 0; i < size; i++) {
    reverse[original[i]] = i;
  }
}
