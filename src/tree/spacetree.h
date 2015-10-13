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
 * @file spacetree.h
 *
 * Generalized space partitioning tree.
 *
 * @experimental
 */

#ifndef TREE_SPACETREE_H
#define TREE_SPACETREE_H

#include "../base/base.h"
#include "statistic.h"
//#include "statistic.h"

/**
 * A binary space partitioning tree, such as KD or ball tree.
 *
 * This particular tree forbids you from having more children.
 *
 * @param TBound the bounding type of each child (TODO explain interface)
 * @param TDataset the data set type
 * @param TStatistic extra data in the node
 *
 * @experimental
 */
template<class TBound,
         class TDataset,
         class TStatistic = EmptyStatistic<TDataset> >
class BinarySpaceTree {
 public:
  typedef TBound Bound;
  typedef TDataset Dataset;
  typedef TStatistic Statistic;

 private:
  Bound bound_;
  BinarySpaceTree *left_;
  BinarySpaceTree *right_;
  index_t begin_;
  index_t count_;
  Statistic stat_;

  OT_DEF(BinarySpaceTree) {
    OT_MY_OBJECT(bound_);
    OT_PTR_NULLABLE(left_);
    OT_PTR_NULLABLE(right_);
    OT_MY_OBJECT(begin_);
    OT_MY_OBJECT(count_);
    OT_MY_OBJECT(stat_);
  }

 public:
  /*
  BinarySpaceTree() {
    DEBUG_ONLY(begin_ = BIG_BAD_NUMBER);
    DEBUG_ONLY(count_ = BIG_BAD_NUMBER);
    DEBUG_POISON_PTR(left_);
    DEBUG_POISON_PTR(right_);
  }

  ~BinarySpaceTree() {
    if (!is_leaf()) {
      delete left_;
      delete right_;
    }
    DEBUG_ONLY(begin_ = BIG_BAD_NUMBER);
    DEBUG_ONLY(count_ = BIG_BAD_NUMBER);
    DEBUG_POISON_PTR(left_);
    DEBUG_POISON_PTR(right_);
  }
  */

  void Init(index_t begin_in, index_t count_in) {
    //DEBUG_ASSERT(begin_ == BIG_BAD_NUMBER);
    DEBUG_POISON_PTR(left_);
    DEBUG_POISON_PTR(right_);
    begin_ = begin_in;
    count_ = count_in;
  }

  /**
   * Find a node in this tree by its begin and count.
   *
   * Every node is uniquely identified by these two numbers.
   * This is useful for communicating position over the network,
   * when pointers would be invalid.
   *
   * @param begin_q the begin() of the node to find
   * @param count_q the count() of the node to find
   * @return the found node, or NULL
   */
  const BinarySpaceTree* FindByBeginCount(
      index_t begin_q, index_t count_q) const {
    DEBUG_ASSERT(begin_q >= begin_);
    DEBUG_ASSERT(count_q <= count_);
    if (begin_ == begin_q && count_ == count_q) {
      return this;
    } else if (unlikely(is_leaf())) {
      return NULL;
    } else if (begin_q < right_->begin_) {
      return left_->FindByBeginCount(begin_q, count_q);
    } else {
      return right_->FindByBeginCount(begin_q, count_q);
    }
  }
  
  /**
   * Find a node in this tree by its begin and count (const).
   *
   * Every node is uniquely identified by these two numbers.
   * This is useful for communicating position over the network,
   * when pointers would be invalid.
   *
   * @param begin_q the begin() of the node to find
   * @param count_q the count() of the node to find
   * @return the found node, or NULL
   */
  BinarySpaceTree* FindByBeginCount(
      index_t begin_q, index_t count_q) {
    DEBUG_ASSERT(begin_q >= begin_);
    DEBUG_ASSERT(count_q <= count_);
    if (begin_ == begin_q && count_ == count_q) {
      return this;
    } else if (unlikely(is_leaf())) {
      return NULL;
    } else if (begin_q < right_->begin_) {
      return left_->FindByBeginCount(begin_q, count_q);
    } else {
      return right_->FindByBeginCount(begin_q, count_q);
    }
  }
  
  // TODO: Not const correct
  
  /**
   * Used only when constructing the tree.
   */
  void set_children(const Dataset& data,
      BinarySpaceTree *left_in, BinarySpaceTree *right_in) {
    left_ = left_in;
    right_ = right_in;
    if (!is_leaf()) {
      stat_.Init(data, begin_, count_, left_->stat_, right_->stat_);
      DEBUG_ASSERT(count_ == left_->count_ + right_->count_);
      DEBUG_ASSERT(left_->begin_ == begin_);
      DEBUG_ASSERT(right_->begin_ == begin_ + left_->count_);
    } else {
      stat_.Init(data, begin_, count_);
    }
  }

  const Bound& bound() const {
    return bound_;
  }

  Bound& bound() {
    return bound_;
  }

  const Statistic& stat() const {
    return stat_;
  }

  Statistic& stat() {
    return stat_;
  }

  bool is_leaf() const {
    return !left_;
  }

  /**
   * Gets the left branch of the tree.
   */
  BinarySpaceTree *left() const {
    // TODO: Const correctness
    return left_;
  }

  /**
   * Gets the right branch.
   */
  BinarySpaceTree *right() const {
    // TODO: Const correctness
    return right_;
  }

  /**
   * Gets the index of the begin point of this subset.
   */
  index_t begin() const {
    return begin_;
  }

  /**
   * Gets the index one beyond the last index in the series.
   */
  index_t end() const {
    return begin_ + count_;
  }
  
  /**
   * Gets the number of points in this subset.
   */
  index_t count() const {
    return count_;
  }
  
  void Print() const {
    printf("node: %d to %d: %d points total\n",
       begin_, begin_ + count_ - 1, count_);
    if (!is_leaf()) {
      left_->Print();
      right_->Print();
    }
  }

};

#endif
