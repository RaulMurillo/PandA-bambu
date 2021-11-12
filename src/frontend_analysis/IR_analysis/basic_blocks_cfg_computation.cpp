/*
 *
 *                   _/_/_/    _/_/   _/    _/ _/_/_/    _/_/
 *                  _/   _/ _/    _/ _/_/  _/ _/   _/ _/    _/
 *                 _/_/_/  _/_/_/_/ _/  _/_/ _/   _/ _/_/_/_/
 *                _/      _/    _/ _/    _/ _/   _/ _/    _/
 *               _/      _/    _/ _/    _/ _/_/_/  _/    _/
 *
 *             ***********************************************
 *                              PandA Project
 *                     URL: http://panda.dei.polimi.it
 *                       Politecnico di Milano - DEIB
 *                        System Architectures Group
 *             ***********************************************
 *              Copyright (C) 2004-2021 Politecnico di Milano
 *
 *   This file is part of the PandA framework.
 *
 *   The PandA framework is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
/**
 * @file basic_blocks_cfg_computation.cpp
 * @brief Build basic block control flow graph data structure starting from the tree_manager.
 *
 * @author Marco Lattuada <lattuada@elet.polimi.it>
 *
 */

/// Autoheader include
#include "config_HAVE_BAMBU_BUILT.hpp"
#include "config_HAVE_ZEBU_BUILT.hpp"

/// Header include
#include "basic_blocks_cfg_computation.hpp"

///. include
#include "Parameter.hpp"

/// behavior includes
#include "application_manager.hpp"
#include "basic_block.hpp"
#include "basic_blocks_graph_constructor.hpp"
#include "function_behavior.hpp"
#include "op_graph.hpp"
#if HAVE_HOST_PROFILING_BUILT
#include "profiling_information.hpp"
#endif
#include "operations_graph_constructor.hpp"

/// design_flow includes
#include "design_flow_graph.hpp"
#include "design_flow_manager.hpp"

/// frontend_analysis/IR_analysis
#include "operations_cfg_computation.hpp"

/// graph include
#include "graph.hpp"

/// STD include
#include <fstream>

/// STL include
#include <list>

/// parser/compiler include
#include "token_interface.hpp"

/// tree include
#include "behavioral_helper.hpp"
#include "behavioral_writer_helper.hpp"
#include "ext_tree_node.hpp"
#include "tree_basic_block.hpp"
#include "tree_manager.hpp"
#include "tree_reindex.hpp"

/// utility include
#include "dbgPrintHelper.hpp"
#include "string_manipulation.hpp" // for GET_CLASS

BasicBlocksCfgComputation::BasicBlocksCfgComputation(const ParameterConstRef _parameters, const application_managerRef _AppM, unsigned int _function_id, const DesignFlowManagerConstRef _design_flow_manager)
    : FunctionFrontendFlowStep(_AppM, _function_id, BASIC_BLOCKS_CFG_COMPUTATION, _design_flow_manager, _parameters)
{
   debug_level = parameters->get_class_debug_level(GET_CLASS(*this), DEBUG_LEVEL_NONE);
}

BasicBlocksCfgComputation::~BasicBlocksCfgComputation() = default;

const CustomUnorderedSet<std::pair<FrontendFlowStepType, FrontendFlowStep::FunctionRelationship>> BasicBlocksCfgComputation::ComputeFrontendRelationships(const DesignFlowStep::RelationshipType relationship_type) const
{
   CustomUnorderedSet<std::pair<FrontendFlowStepType, FunctionRelationship>> relationships;
   switch(relationship_type)
   {
      case(DEPENDENCE_RELATIONSHIP):
      {
         relationships.insert(std::make_pair(COMPLETE_BB_GRAPH, SAME_FUNCTION));
#if HAVE_BAMBU_BUILT
         if(parameters->isOption(OPT_writer_language))
         {
            relationships.insert(std::make_pair(HDL_VAR_DECL_FIX, SAME_FUNCTION));
         }
         else
#endif
         {
            relationships.insert(std::make_pair(VAR_DECL_FIX, SAME_FUNCTION));
         }
#if HAVE_BAMBU_BUILT
         relationships.insert(std::make_pair(EXTRACT_PATTERNS, SAME_FUNCTION));
#endif
         break;
      }
      case(INVALIDATION_RELATIONSHIP):
      {
         break;
      }
      case(PRECEDENCE_RELATIONSHIP):
      {
#if HAVE_ZEBU_BUILT
         relationships.insert(std::make_pair(SHORT_CIRCUIT_STRUCTURING, SAME_FUNCTION));
         relationships.insert(std::make_pair(SPLIT_PHINODES, SAME_FUNCTION));
#endif
#if HAVE_BAMBU_BUILT
         relationships.insert(std::make_pair(COND_EXPR_RESTRUCTURING, SAME_FUNCTION));
         relationships.insert(std::make_pair(CSE_STEP, SAME_FUNCTION));
         relationships.insert(std::make_pair(FANOUT_OPT, SAME_FUNCTION));
         relationships.insert(std::make_pair(FUNCTION_CALL_OPT, SAME_FUNCTION));
         relationships.insert(std::make_pair(IR_LOWERING, SAME_FUNCTION));
#if HAVE_ILP_BUILT && HAVE_BAMBU_BUILT
         relationships.insert(std::make_pair(SDC_CODE_MOTION, SAME_FUNCTION));
#endif
         relationships.insert(std::make_pair(UN_COMPARISON_LOWERING, SAME_FUNCTION));
#endif
         break;
      }
      default:
      {
         THROW_UNREACHABLE("");
      }
   }
   return relationships;
}

void BasicBlocksCfgComputation::Initialize()
{
   if(bb_version != 0 and bb_version != function_behavior->GetBBVersion())
   {
      function_behavior->bbgc->Clear();
#if HAVE_HOST_PROFILING_BUILT
      function_behavior->profiling_information->Clear();
#endif
   }
}

DesignFlowStep_Status BasicBlocksCfgComputation::InternalExec()
{
   const tree_managerRef TM = AppM->get_tree_manager();
   const BasicBlocksGraphConstructorRef bbgc = function_behavior->bbgc;
   tree_nodeRef tn = TM->get_tree_node_const(function_id);
   auto* fd = GetPointer<function_decl>(tn);
   THROW_ASSERT(fd && fd->body, "Node is not a function or it hasn't a body");
   auto* sl = GetPointer<statement_list>(GET_NODE(fd->body));
   THROW_ASSERT(sl, "Body is not a statement_list");
   std::map<unsigned int, blocRef>::iterator it_bb, it_bb_end = sl->list_of_bloc.end();
   for(it_bb = sl->list_of_bloc.begin(); it_bb != it_bb_end; ++it_bb)
   {
      if(it_bb->second->number != BB_ENTRY and it_bb->second->number != BB_EXIT)
      {
         continue;
      }
      bbgc->add_vertex(it_bb->second);
      INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "---Added basic block with index " + boost::lexical_cast<std::string>(it_bb->second->number));
      if(it_bb->second->number == BB_EXIT)
      {
         const vertex exit = bbgc->Cget_vertex(BB_EXIT);
         bbgc->connect_to_entry(exit);
      }
   }
   for(it_bb = sl->list_of_bloc.begin(); it_bb != it_bb_end; ++it_bb)
   {
      if(it_bb->second->number == BB_ENTRY || it_bb->second->number == BB_EXIT)
      {
         continue;
      }
      bbgc->add_vertex(it_bb->second);
      INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "---Added basic block with index " + boost::lexical_cast<std::string>(it_bb->second->number));
   }
   auto b_end = sl->list_of_bloc.end();
   for(auto b = sl->list_of_bloc.begin(); b != b_end; ++b)
   {
      if(b->second->number == BB_ENTRY || b->second->number == BB_EXIT)
      {
         continue;
      }
      INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "-->Considering connections for BB" + boost::lexical_cast<std::string>(b->first));
      const vertex current = bbgc->Cget_vertex(b->second->number);
      if(b->second->list_of_pred.empty() or std::find(b->second->list_of_pred.begin(), b->second->list_of_pred.end(), bloc::ENTRY_BLOCK_ID) != b->second->list_of_pred.end())
      {
         INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "---Connecting Basic block " + boost::lexical_cast<std::string>(b->second->number) + " to ENTRY");
         bbgc->connect_to_entry(current);
      }
      if(b->second->list_of_succ.empty())
      {
         // PRINT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "Connecting Basic block " + boost::lexical_cast<std::string>((*b)->number) + " to EXIT");
         // bbgc->connect_to_exit(current);
      }
      else
      {
         auto su_end = b->second->list_of_succ.end();
         for(auto su = b->second->list_of_succ.begin(); su != su_end; ++su)
         {
            if((*su) == bloc::EXIT_BLOCK_ID)
            {
               INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "---Connecting Basic block " + boost::lexical_cast<std::string>(b->second->number) + " to EXIT");
               bbgc->connect_to_exit(current);
            }
            else
            {
               INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "---Connecting Basic block " + boost::lexical_cast<std::string>(b->second->number) + " to " + boost::lexical_cast<std::string>(*su));
               bbgc->AddEdge(current, bbgc->Cget_vertex(*su), CFG_SELECTOR);
               // Considering label
               if(*su == b->second->true_edge)
               {
                  bbgc->add_bb_edge_info(current, bbgc->Cget_vertex(*su), CFG_SELECTOR, T_COND);
               }
               if(*su == b->second->false_edge)
               {
                  bbgc->add_bb_edge_info(current, bbgc->Cget_vertex(*su), CFG_SELECTOR, F_COND);
               }
            }
         }
         const auto& statements = b->second->CGetStmtList();
         if(!statements.empty())
         {
            tree_nodeRef last = statements.back();
            /// switch statements
            if(GET_NODE(last)->get_kind() == gimple_switch_K)
            {
               INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "-->BB" + boost::lexical_cast<std::string>(b->first) + " ends with a switch");
               // Map between gimple_label and index of basic block
               std::map<tree_nodeRef, unsigned int> label_to_bb;
               su_end = b->second->list_of_succ.end();
               for(auto su = b->second->list_of_succ.begin(); su != su_end; ++su)
               {
                  THROW_ASSERT(sl->list_of_bloc[*su]->CGetStmtList().size(), "Empty Basic Block");
                  const auto first = sl->list_of_bloc[*su]->CGetStmtList().front();
                  THROW_ASSERT(GetPointer<gimple_label>(GET_NODE(first)), "First operation of BB" + STR(*su) + " is a " + GET_NODE(first)->get_kind_text() + ": " + GET_NODE(first)->ToString());
                  auto* le = GetPointer<gimple_label>(GET_NODE(first));
                  label_to_bb[GET_NODE(le->op)] = *su;
                  INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "---Gimple label of BB" + boost::lexical_cast<std::string>(*su) + " is " + boost::lexical_cast<std::string>(GET_INDEX_NODE(le->op)));
               }
               auto* se = GetPointer<gimple_switch>(GET_NODE(last));
               THROW_ASSERT(se->op1, "case_label_exprs not found");
               auto* tv = GetPointer<tree_vec>(GET_NODE(se->op1));
               auto it_end = tv->list_of_op.end();
               for(auto it = tv->list_of_op.begin(); it != it_end; ++it)
               {
                  auto* cl = GetPointer<case_label_expr>(GET_NODE(*it));
                  THROW_ASSERT(label_to_bb.find(GET_NODE(cl->got)) != label_to_bb.end(), "There is not corresponding case_label_exprs with index " + boost::lexical_cast<std::string>(GET_INDEX_NODE(cl->got)));
                  if(cl->default_flag)
                  {
                     bbgc->add_bb_edge_info(current, bbgc->Cget_vertex(label_to_bb[GET_NODE(cl->got)]), CFG_SELECTOR, default_COND);
                  }
                  else
                  {
                     bbgc->add_bb_edge_info(current, bbgc->Cget_vertex(label_to_bb[GET_NODE(cl->got)]), CFG_SELECTOR, GET_INDEX_NODE(*it));
                  }
               }
               INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "<--");
            }
            /// computed goto
            else if(GET_NODE(last)->get_kind() == gimple_goto_K && b->second->list_of_succ.size() > 1)
            {
               // Map between gimple_label and index of basic block
               su_end = b->second->list_of_succ.end();
               for(auto su = b->second->list_of_succ.begin(); su != su_end; ++su)
               {
                  bbgc->add_bb_edge_info(current, bbgc->Cget_vertex(*su), CFG_SELECTOR, *su);
               }
            }
            /// multi-way if
            else if(GET_NODE(last)->get_kind() == gimple_multi_way_if_K)
            {
               auto* gmwi = GetPointer<gimple_multi_way_if>(GET_NODE(last));
               for(const auto& cond : gmwi->list_of_cond)
               {
                  bbgc->add_bb_edge_info(current, bbgc->Cget_vertex(cond.second), CFG_SELECTOR, cond.first ? cond.first->index : default_COND);
               }
            }
         }
      }
      INDENT_DBG_MEX(DEBUG_LEVEL_VERY_PEDANTIC, debug_level, "<--Considered connections for BB" + boost::lexical_cast<std::string>(b->first));
   }
   const vertex exit = bbgc->Cget_vertex(BB_EXIT);
   const BBGraphRef fcfg = function_behavior->GetBBGraph(FunctionBehavior::FBB);
   VertexIterator v, v_end;
   for(boost::tie(v, v_end) = boost::vertices(*fcfg); v != v_end; v++)
   {
      if(boost::out_degree(*v, *fcfg) == 0 and *v != exit)
      {
         bbgc->AddEdge(*v, exit, CFG_SELECTOR);
      }
   }
   if(parameters->getOption<bool>(OPT_print_dot))
   {
      function_behavior->GetBBGraph(FunctionBehavior::BB)->WriteDot("BB_CFG.dot");
   }
   return DesignFlowStep_Status::SUCCESS;
}
