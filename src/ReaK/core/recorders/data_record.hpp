/**
 * \file data_record.hpp
 *
 * This library declares the basic class for data recording to a file. Here, "data" is meant as
 * columns of floating-point records of data, such as simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date january 2010
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef DATA_RECORD_HPP
#define DATA_RECORD_HPP

#include "base/defs.hpp"

#include "base/thread_incl.hpp"



#include <string>
#include <exception>
#include <vector>
#include <queue>
#include <map>
#include <iostream>

#include "base/shared_object.hpp"

#include "rtti/so_type.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK's Data Recorders and Extractors */
namespace recorder {

/**
 * This exception is thrown whenever a data entry is written passed the
 * a-priori fixed column count.
 */
class out_of_bounds : public std::exception {
  public:

    out_of_bounds() { };

    ~out_of_bounds() throw() {};

    const char* what() const throw() {
      return "Data record went out of bounds!";
    };
};

/**
 * This exception is thrown whenever a data entry is read passed the
 * number of rows in the record.
 */
class end_of_record : public std::exception {
  public:

    end_of_record() { };

    ~end_of_record() throw() {};

    const char* what() const throw() {
      return "No more data rows in the record!";
    };
};

/**
 * This exception is thrown whenever a special flag was given to the data
 * recorder that is not valid, such as terminating the column name list
 * when it has already been closed.
 */
class improper_flag : public std::exception {
  public:

    improper_flag() { };

    ~improper_flag() throw() {};

    const char* what() const throw() {
      return "Flagged operation requested on data recorder was invalid!";
    };
};


/* forward-declarations */
class data_recorder;
class data_extractor;



/**
 * This class is used to represent a complete row of entries (values) for a data recorder or extractor, 
 * each associated with a column-name (names must be unique).
 * \note Objects of this class are created using the 'getFreshNamedValueRow' function in the data recorder / extractor classes.
 */
class named_value_row {
  private:
    const std::map<std::string, std::size_t>* p_named_indices;
    std::vector<double> values;
    
    explicit named_value_row(const std::map<std::string, std::size_t>& aMap) : p_named_indices(&aMap), values(aMap.size(), 0.0) { };
    
  public:
    
    friend class data_recorder;
    friend class data_extractor;
    
    /**
     * Entry-access function. This function can be used to obtain write-access to an entry associated to a given name.
     * \param s The name of the entry addressed.
     * \return A reference to the entry corresponding to the given name.
     */
    double& operator[](const std::string& s) {
      std::map<std::string, std::size_t>::const_iterator it = p_named_indices->find(s);
      if( it == p_named_indices->end() )
        throw out_of_bounds();
      return values[ it->second ];
    };
    /**
     * Entry-access function. This function can be used to obtain read-access to an entry associated to a given name.
     * \param s The name of the entry addressed.
     * \return The entry corresponding to the given name.
     */
    double operator[](const std::string& s) const {
      std::map<std::string, std::size_t>::const_iterator it = p_named_indices->find(s);
      if( it == p_named_indices->end() )
        throw out_of_bounds();
      return values[ it->second ];
    };
    
};
  



/**
 * This class is the basis for all data recording classes. This class handles the basic
 * operations for buffering of the data and column name records.
 */
class data_recorder : public shared_object {
  protected:
    volatile unsigned int colCount; ///< Holds the column count.
    volatile unsigned int rowCount; ///< Holds the number of rows of data records.
    volatile unsigned int currentColumn; ///< Holds the current column to which the next data entry will be written to.
    unsigned int flushSampleRate; ///< Holds the sample rate at which the data is automatically flushed to the file.
    unsigned int maxBufferSize; ///< Holds the maximum size for the data buffer, overload will trigger a file-flush.
    std::vector<std::string> names; ///< Holds the list of column names.
    mutable std::map<std::string, std::size_t> named_indices; ///< Holds the map from the column names to the index within a value-row.
    std::queue<double> values_rm; ///< Holds the data buffer.
    shared_ptr<std::ostream> out_stream; ///< Holds the output-stream of the data record.
    
    ReaKaux::mutex access_mutex; ///< Mutex to lock the read/write on the data buffer.
    ReaK::shared_ptr<ReaKaux::thread> writing_thread; ///< Holds the instance of the data writing thread.
    
    /**
     * This class is used as a callable function-object for data writing thread.
     */
    struct record_process {
      public:
        data_recorder& parent;
        record_process(	data_recorder& aParent ) : parent(aParent) { };
        void operator()();
    };
    
    /**
     * Overridable function which writes data to the file in whichever format specific to the derived class.
     */
    virtual void writeRow() { };
    /**
     * Overridable function which writes column names to the file in whichever format specific to the derived class.
     */
    virtual void writeNames() { };
    
    virtual void setStreamImpl(const shared_ptr<std::ostream>& aStreamPtr) = 0;
    
  public:
    
    /**
     * This function is the factory to create named-value-row objects to represent a row of entries, addressable by name.
     * \return A fresh object that is ready to accept all the values of a row of entries to the data recorder.
     */
    named_value_row getFreshNamedValueRow() const {
      return named_value_row(named_indices);
    };
    
    /**
     * This function returns the number of columns in this recorder.
     * \return The number of columns in this recorder.
     */
    unsigned int getColCount() const { return colCount; };
    
    /// Data record-specific flags for special operations.
    enum flag {
      end_name_row, ///< Ends the definition of the column name list.
      end_value_row, ///< Ends the recording of the data entry for the current row.
      flush, ///< Flushes the data buffer to the file.
      close ///< Closes the file and data buffer is flushed.
    };
    
    /**
     * Default Constructor.
     */
    data_recorder() : shared_object(),
                      colCount(0),
                      rowCount(0),
                      currentColumn(0),
                      flushSampleRate(50),
                      maxBufferSize(500),
                      names(),
                      values_rm(),
                      out_stream(),
                      access_mutex(),
                      writing_thread() { };
    
    /**
     * Destructor.
     */
    virtual ~data_recorder();
    
    /**
     * Operator to record a data entry.
     */
    data_recorder& operator <<(double value);
    
    /**
     * Operator to record an entire row of named entries.
     */
    data_recorder& operator <<(const named_value_row& values);
    
    /**
     * Operator to record a column name.
     */
    data_recorder& operator <<(const std::string& name);
    
    /**
     * Operator to record a column name.
     */
    data_recorder& operator <<(const char* name) { return (*this) << std::string(name); };
    
    /**
     * Operator to issue a special operation's flag.
     */
    data_recorder& operator <<(flag some_flag);
    
    /**
     * Operator to record a vector of in-order (nameless) entries.
     */
    template <typename Vector>
    data_recorder& operator <<(const Vector& values) {
      for(std::size_t i = 0; i < values.size(); ++i)
        (*this) << values[i];
      return *this;
    };
    
    /**
     * Sets the stream.
     */
    void setStream(std::ostream& aStream) {
      setStreamImpl(shared_ptr<std::ostream>(&aStream, null_deleter()));
    };
    
    /**
     * Sets the stream via a shared-pointer.
     */
    void setStream(const shared_ptr<std::ostream>& aStreamPtr) {
      setStreamImpl(aStreamPtr);
    };
    
    /**
     * Sets the filename for the file.
     */
    virtual void setFileName(const std::string& aFilename);
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(static_cast<unsigned int>(colCount))
        & RK_SERIAL_SAVE_WITH_NAME(flushSampleRate)
        & RK_SERIAL_SAVE_WITH_NAME(maxBufferSize)
        & RK_SERIAL_SAVE_WITH_NAME(names);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
      ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
      colCount = 0;
      if(writing_thread) {
        lock_here.unlock();
        writing_thread->join();
        lock_here.lock();
      };
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      unsigned int aColCount;
      A & RK_SERIAL_LOAD_WITH_ALIAS("colCount",aColCount)
        & RK_SERIAL_LOAD_WITH_NAME(flushSampleRate)
        & RK_SERIAL_LOAD_WITH_NAME(maxBufferSize)
        & RK_SERIAL_LOAD_WITH_NAME(names);
      colCount = aColCount;
      for(std::size_t i = 0; i < colCount; ++i)
        named_indices[names[i]] = i;
      rowCount = 0;
      currentColumn = 0;
      values_rm = std::queue<double>();
      lock_here.unlock();
      writing_thread = ReaK::shared_ptr<ReaKaux::thread>(new ReaKaux::thread(record_process(*this)));
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(data_recorder,0x81100001,1,"data_recorder",shared_object)
    
};



/**
 * This class is the basis for all data extracting classes. This class handles the basic
 * operations for buffering of the data and column name records.
 */
class data_extractor : public shared_object {
  protected:
    volatile unsigned int colCount; ///< Holds the column count.
    volatile unsigned int currentColumn; ///< Holds the current column to which the next data entry will be read from.
    volatile unsigned int currentNameCol; ///< Holds the current column to which the next name entry will be read from.
    unsigned int flushSampleRate; ///< Holds the sample rate at which the data is automatically flushed to the file.
    unsigned int minBufferSize; ///< Holds the minimum size for the data buffer, underload will trigger a file-read.
    std::vector<std::string> names; ///< Holds the list of column names.
    mutable std::map<std::string, std::size_t> named_indices; ///< Holds the map from the column names to the index within a value-row.
    std::queue<double> values_rm; ///< Holds the data buffer.
    shared_ptr<std::istream> in_stream; ///< Holds the input-stream of the data record.
    
    ReaKaux::mutex access_mutex; ///< Mutex to lock the read/write on the data buffer.
    ReaK::shared_ptr<ReaKaux::thread> reading_thread; ///< Holds the instance of the data writing thread.

    /**
     * This class is used as a callable function-object for data writing thread.
     */
    struct extract_process {
      public:
        data_extractor& parent;
        extract_process(data_extractor& aParent ) : parent(aParent) { };
        void operator()();
    };

    /**
     * Overridable function which writes data to the file in whichever format specific to the derived class.
     */
    virtual bool readRow() { return true; };
    /**
     * Overridable function which writes column names to the file in whichever format specific to the derived class.
     */
    virtual bool readNames() { return true; };
    /**
     * Sets the filename for the file, overridable for the derived class which handles the file IO (or other).
     */
    virtual void setStreamImpl(const shared_ptr<std::istream>& aStreamPtr) = 0;
    
  public:
    
    /**
     * This function is the factory to create named-value-row objects to represent a row of entries, addressable by name.
     * \return A fresh object that is ready to accept all the values of a row of entries to the data recorder.
     */
    named_value_row getFreshNamedValueRow() const {
      if( named_indices.size() != names.size() ) {
        for(std::size_t i = 0; i < names.size(); ++i)
          named_indices[names[i]] = i;
      };
      return named_value_row(named_indices);
    };
    
    /**
     * This function returns the number of columns in this extractor.
     * \return The number of columns in this extractor.
     */
    unsigned int getColCount() const { return colCount; };

    /// Data record-specific flags for special operations.
    enum flag {
      end_value_row, ///< Ends the recording of the data entry for the current row.
      advance, ///< Reads the data buffer from the file.
      close ///< Closes the file and data buffer is flushed.
    };

    /**
     * Default Constructor.
     */
    data_extractor() : shared_object(),
                      colCount(0),
                      currentColumn(0),
                      currentNameCol(0),
                      flushSampleRate(50),
                      minBufferSize(20),
                      names(),
                      values_rm(),
                      in_stream(),
                      access_mutex(),
                      reading_thread() { };
    
    /**
     * Destructor.
     */
    virtual ~data_extractor();
    
    /**
     * Operator to extract a data entry.
     */
    data_extractor& operator >>(double& value);
    
    /**
     * Operator to extract a column name.
     */
    data_extractor& operator >>(std::string& name);
    
    /**
     * Operator to extract an entire row of named entries.
     */
    data_extractor& operator >>(named_value_row& values);
    
    /**
     * Operator to issue a special operation's flag.
     */
    data_extractor& operator >>(flag some_flag);
    
    /**
     * Operator to extract a vector of in-order (nameless) entries.
     */
    template <typename Vector>
    data_extractor& operator >>(Vector& values) {
      for(std::size_t i = 0; i < values.size(); ++i)
        (*this) >> values[i];
      return *this;
    };
    
    /**
     * Sets the stream.
     */
    void setStream(std::istream& aStream) {
      setStreamImpl(shared_ptr<std::istream>(&aStream, null_deleter()));
    };
    
    /**
     * Sets the stream via a shared-pointer.
     */
    void setStream(const shared_ptr<std::istream>& aStreamPtr) {
      setStreamImpl(aStreamPtr);
    };
    
    /**
     * Sets the filename for the file, overridable for the derived class which handles the file IO (or other).
     */
    virtual void setFileName(const std::string& aFileName);
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(static_cast<unsigned int>(colCount))
        & RK_SERIAL_SAVE_WITH_NAME(flushSampleRate)
        & RK_SERIAL_SAVE_WITH_NAME(minBufferSize)
        & RK_SERIAL_SAVE_WITH_NAME(names);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
      colCount = 0;
      if(reading_thread) {
        lock_here.unlock();
        reading_thread->join();
        lock_here.lock();
      };
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      unsigned int aColCount;
      A & RK_SERIAL_LOAD_WITH_ALIAS("colCount",aColCount)
        & RK_SERIAL_LOAD_WITH_NAME(flushSampleRate)
        & RK_SERIAL_LOAD_WITH_NAME(minBufferSize)
        & RK_SERIAL_LOAD_WITH_NAME(names);
      colCount = aColCount;
      for(std::size_t i = 0; i < names.size(); ++i)
        named_indices[names[i]] = i;
      currentColumn = 0;
      currentNameCol = 0;
      values_rm = std::queue<double>();
      lock_here.unlock();
      reading_thread = ReaK::shared_ptr<ReaKaux::thread>(new ReaKaux::thread(extract_process(*this)));
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(data_extractor,0x81200001,1,"data_extractor",shared_object)

};


};


};


#endif














