/**
 * \file objtree_archiver.hpp
 *
 * This library declares the class for a creating type schemes that represent the fields contained in the 
 * serialization of a given type.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date November 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
 *
 *    THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
 *
 *    This file is part of ReaK.
 *
 *    ReaK is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ReaK is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ReaK (as LICENSE in the root folder).  
 *    If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef REAK_OBJTREE_ARCHIVER_HPP
#define REAK_OBJTREE_ARCHIVER_HPP

#include "archiver.hpp"

#include "base/serializable.hpp"

#include <fstream>

#include <boost/graph/adjacency_list.hpp>

namespace ReaK {

namespace serialization {

struct object_graph_node {
  shared_ptr< serializable > p_obj;
  std::string xml_src;
};

typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, object_graph_node > object_graph;


/**
 * An input archive used to generate an object tree.
 */
class objtree_iarchive : public iarchive {
  private:
    shared_ptr< object_graph > obj_graph;
    shared_ptr< std::stringstream > current_ss;

    char getNextChar();
    std::string readToken();
    void skipToEndToken(const std::string& name);
    void trimStr(std::string& s);
    bool readNamedValue(const std::string& value_name,std::string& value_str);
    archive_object_header readHeader(const std::string& obj_name);

  protected:

    virtual iarchive& RK_CALL load_serializable_ptr(serializable_shared_pointer& Item);

    virtual iarchive& RK_CALL load_serializable_ptr(const std::pair<std::string, serializable_shared_pointer& >& Item);

    virtual iarchive& RK_CALL load_serializable(serializable& Item);

    virtual iarchive& RK_CALL load_serializable(const std::pair<std::string, serializable& >& Item);

    virtual iarchive& RK_CALL load_char(char& i);

    virtual iarchive& RK_CALL load_char(const std::pair<std::string, char& >& i);

    virtual iarchive& RK_CALL load_unsigned_char(unsigned char& u);

    virtual iarchive& RK_CALL load_unsigned_char(const std::pair<std::string, unsigned char& >& u);

    virtual iarchive& RK_CALL load_int(int& i);

    virtual iarchive& RK_CALL load_int(const std::pair<std::string, int& >& i);

    virtual iarchive& RK_CALL load_unsigned_int(unsigned int& u);

    virtual iarchive& RK_CALL load_unsigned_int(const std::pair<std::string, unsigned int& >& u);

    virtual iarchive& RK_CALL load_float(float& f);

    virtual iarchive& RK_CALL load_float(const std::pair<std::string, float& >& f);

    virtual iarchive& RK_CALL load_double(double& d);

    virtual iarchive& RK_CALL load_double(const std::pair<std::string, double& >& d);

    virtual iarchive& RK_CALL load_bool(bool& b);

    virtual iarchive& RK_CALL load_bool(const std::pair<std::string, bool& >& b);

    virtual iarchive& RK_CALL load_string(std::string& s);

    virtual iarchive& RK_CALL load_string(const std::pair<std::string, std::string& >& s);
  
  public:

    objtree_iarchive(const shared_ptr< object_graph >& aObjGraph);
    virtual ~objtree_iarchive();

};


/**
 * XML output archive.
 */
class objtree_oarchive : public oarchive {
  private:
    shared_ptr< object_graph > obj_graph;
    shared_ptr< std::stringstream > current_ss;
    boost::graph_traits< object_graph >::vertex_descriptor current_node;
    
  protected:
    
    virtual oarchive& RK_CALL saveToNewArchive_impl(const serializable_shared_pointer& Item, const std::string& FileName);
    
    virtual oarchive& RK_CALL saveToNewArchiveNamed_impl(const std::pair<std::string, const serializable_shared_pointer& >& Item, const std::string& FileName);
    
    virtual oarchive& RK_CALL save_serializable_ptr(const serializable_shared_pointer& Item);
    
    virtual oarchive& RK_CALL save_serializable_ptr(const std::pair<std::string, const serializable_shared_pointer& >& Item);
    
    virtual oarchive& RK_CALL save_serializable(const serializable& Item);
    
    virtual oarchive& RK_CALL save_serializable(const std::pair<std::string, const serializable& >& Item);
    
    virtual oarchive& RK_CALL save_char(char i);
    
    virtual oarchive& RK_CALL save_char(const std::pair<std::string, char >& i);
    
    virtual oarchive& RK_CALL save_unsigned_char(unsigned char u);
    
    virtual oarchive& RK_CALL save_unsigned_char(const std::pair<std::string, unsigned char >& u);
    
    virtual oarchive& RK_CALL save_int(int i);
    
    virtual oarchive& RK_CALL save_int(const std::pair<std::string, int >& i);
    
    virtual oarchive& RK_CALL save_unsigned_int(unsigned int u);
    
    virtual oarchive& RK_CALL save_unsigned_int(const std::pair<std::string, unsigned int >& u);
    
    virtual oarchive& RK_CALL save_float(float f);
    
    virtual oarchive& RK_CALL save_float(const std::pair<std::string, float >& f);
    
    virtual oarchive& RK_CALL save_double(double d);
    
    virtual oarchive& RK_CALL save_double(const std::pair<std::string, double >& d);
    
    virtual oarchive& RK_CALL save_bool(bool b);
    
    virtual oarchive& RK_CALL save_bool(const std::pair<std::string, bool >& b);
    
    virtual oarchive& RK_CALL save_string(const std::string& s);
    
    virtual oarchive& RK_CALL save_string(const std::pair<std::string, const std::string& >& s);
    
  public:
    
    objtree_oarchive(const shared_ptr< object_graph >& aObjGraph);
    virtual ~objtree_oarchive();
    
};



}; //serialization

}; //ReaK

#endif






