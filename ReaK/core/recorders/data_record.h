/**
 * \file data_record.h
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

#ifndef REAK_CORE_RECORDERS_DATA_RECORD_H_
#define REAK_CORE_RECORDERS_DATA_RECORD_H_

#include "ReaK/core/base/shared_object.h"

#include <atomic>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

/** Main namespace for ReaK */
namespace ReaK::recorder {

/**
 * This exception is thrown whenever a data entry is written passed the
 * a-priori fixed column count.
 */
class out_of_bounds : public std::exception {
 public:
  out_of_bounds() = default;

  ~out_of_bounds() noexcept override = default;

  [[nodiscard]] const char* what() const noexcept override {
    return "Data record went out of bounds!";
  }
};

/**
 * This exception is thrown whenever a data entry is read passed the
 * number of rows in the record.
 */
class end_of_record : public std::exception {
 public:
  end_of_record() = default;

  ~end_of_record() noexcept override = default;

  [[nodiscard]] const char* what() const noexcept override {
    return "No more data rows in the record!";
  }
};

/**
 * This exception is thrown whenever a special flag was given to the data
 * recorder that is not valid, such as terminating the column name list
 * when it has already been closed.
 */
class improper_flag : public std::exception {
 public:
  improper_flag() = default;

  ~improper_flag() noexcept override = default;

  [[nodiscard]] const char* what() const noexcept override {
    return "Flagged operation requested on data recorder was invalid!";
  }
};

/* forward-declarations */
class data_recorder;
class data_extractor;

/**
 * This class is used to represent a complete row of entries (values) for a data recorder or extractor,
 * each associated with a column-name (names must be unique).
 * \note Objects of this class are created using the 'getFreshNamedValueRow' function in the data recorder / extractor
 * classes.
 */
class named_value_row {
 private:
  const std::map<std::string, std::size_t>* p_named_indices;
  std::vector<double> values;

  explicit named_value_row(const std::map<std::string, std::size_t>& aMap)
      : p_named_indices(&aMap), values(aMap.size(), 0.0) {}

 public:
  friend class data_recorder;
  friend class data_extractor;

  /**
   * Entry-access function. This function can be used to obtain write-access to an entry associated to a given name.
   * \param s The name of the entry addressed.
   * \return A reference to the entry corresponding to the given name.
   */
  double& operator[](const std::string& s) {
    auto it = p_named_indices->find(s);
    if (it == p_named_indices->end()) {
      throw out_of_bounds();
    }
    return values[it->second];
  }
  /**
   * Entry-access function. This function can be used to obtain read-access to an entry associated to a given name.
   * \param s The name of the entry addressed.
   * \return The entry corresponding to the given name.
   */
  double operator[](const std::string& s) const {
    auto it = p_named_indices->find(s);
    if (it == p_named_indices->end()) {
      throw out_of_bounds();
    }
    return values[it->second];
  }
};

/**
 * This class is the basis for all data recording classes. This class handles the basic
 * operations for buffering of the data and column name records.
 */
class data_recorder : public shared_object {
 protected:
  // Holds the column count.
  std::atomic<unsigned int> colCount;
  // Holds the number of rows of data records.
  std::atomic<unsigned int> rowCount;
  // Holds the current column to which the next data entry will be written to.
  std::atomic<unsigned int> currentColumn;
  // Holds the sample rate at which the data is automatically flushed to the file.
  unsigned int flushSampleRate{50};
  // Holds the maximum size for the data buffer, overload will trigger a file-flush.
  unsigned int maxBufferSize{500};
  // Holds the list of column names.
  std::vector<std::string> names;
  // Holds the map from the column names to the index within a value-row.
  mutable std::map<std::string, std::size_t> named_indices;
  // Holds the data buffer.
  std::queue<double> values_rm;
  // Holds the output-stream of the data record.
  std::shared_ptr<std::ostream> out_stream;

  // Mutex to lock the read/write on the data buffer.
  std::mutex access_mutex;
  // Holds the instance of the data writing thread.
  std::shared_ptr<std::thread> writing_thread;

  /**
   * This class is used as a callable function-object for data writing thread.
   */
  struct record_process {
   public:
    data_recorder& parent;
    explicit record_process(data_recorder& aParent) : parent(aParent){};
    void operator()();
  };

  void closeRecordProcess();

  /**
   * Overridable function which writes data to the file in whichever format specific to the derived class.
   */
  virtual void writeRow() {}
  /**
   * Overridable function which writes column names to the file in whichever format specific to the derived class.
   */
  virtual void writeNames() {}

  virtual void setStreamImpl(
      const std::shared_ptr<std::ostream>& aStreamPtr) = 0;

 public:
  /**
   * Returns the flushing sample rate of this data streamer in Hz.
   * \note A flushing sample rate of 0 signifies immediate transmission.
   * \return The flushing sample rate of this data streamer in Hz.
   */
  unsigned int getFlushSampleRate() const { return flushSampleRate; }

  /**
   * Sets the flushing sample rate of this data streamer in Hz.
   * \note A flushing sample rate of 0 signifies immediate transmission.
   * \param aFlushSampleRate The flushing sample rate of this data streamer in Hz.
   */
  void setFlushSampleRate(unsigned int aFlushSampleRate) {
    flushSampleRate = aFlushSampleRate;
  }

  /**
   * Returns the maximum size of the data buffer, after which, data must be flushed to the stream.
   * \return The maximum size of the data buffer, after which, data must be flushed to the stream.
   */
  unsigned int getMaxBufferSize() const { return maxBufferSize; }

  /**
   * Sets the maximum size of the data buffer, after which, data must be flushed to the stream.
   * \param aMinBufferSize The maximum size of the data buffer, after which, data must be flushed to the stream.
   */
  void setMaxBufferSize(unsigned int aMaxBufferSize) {
    maxBufferSize = aMaxBufferSize;
  }

  /**
   * This function is the factory to create named-value-row objects to represent a row of entries, addressable by name.
   * \return A fresh object that is ready to accept all the values of a row of entries to the data recorder.
   */
  named_value_row getFreshNamedValueRow() const {
    return named_value_row(named_indices);
  }

  /**
   * This function returns the number of columns in this recorder.
   * \return The number of columns in this recorder.
   */
  unsigned int getColCount() const { return colCount; }

  /// Data record-specific flags for special operations.
  enum flag {
    end_name_row,  ///< Ends the definition of the column name list.
    end_value_row,  ///< Ends the recording of the data entry for the current row.
    flush,          ///< Flushes the data buffer to the file.
    close           ///< Closes the file and data buffer is flushed.
  };

  /**
   * Default Constructor.
   */
  data_recorder()
      : shared_object(), colCount(0), rowCount(0), currentColumn(0) {}

  /**
   * Destructor.
   */
  ~data_recorder() override;

  /**
   * Operator to record a data entry.
   */
  data_recorder& operator<<(double value);

  /**
   * Operator to record an entire row of named entries.
   */
  data_recorder& operator<<(const named_value_row& values);

  /**
   * Operator to record a column name.
   */
  data_recorder& operator<<(const std::string& name);

  /**
   * Operator to record a vector of column names.
   */
  data_recorder& operator<<(const std::vector<std::string>& aNames) {
    for (const auto& name : aNames) {
      (*this) << name;
    }
    return *this;
  }

  /**
   * Operator to record a column name.
   */
  data_recorder& operator<<(const char* name) {
    return (*this) << std::string(name);
  }

  /**
   * Operator to issue a special operation's flag.
   */
  data_recorder& operator<<(flag some_flag);

  /**
   * Operator to record a vector of in-order (nameless) entries.
   */
  template <typename Vector>
  data_recorder& operator<<(const Vector& values) {
    for (std::size_t i = 0; i < values.size(); ++i) {
      (*this) << values[i];
    }
    return *this;
  }

  /**
   * Sets the stream.
   */
  void setStream(std::ostream& aStream) {
    setStreamImpl(std::shared_ptr<std::ostream>(&aStream, null_deleter()));
  }

  /**
   * Sets the stream via a shared-pointer.
   */
  void setStream(const std::shared_ptr<std::ostream>& aStreamPtr) {
    setStreamImpl(aStreamPtr);
  }

  /**
   * Sets the filename for the file.
   */
  virtual void setFileName(const std::string& aFilename);

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override;
  void load(serialization::iarchive& A, unsigned int /*Version*/) override;

  RK_RTTI_MAKE_ABSTRACT_1BASE(data_recorder, 0x81100001, 1, "data_recorder",
                              shared_object)
};

/**
 * This class is the basis for all data extracting classes. This class handles the basic
 * operations for buffering of the data and column name records.
 */
class data_extractor : public shared_object {
 protected:
  // Holds the column count.
  std::atomic<unsigned int> colCount;
  // Holds the current column to which the next data entry will be read from.
  std::atomic<unsigned int> currentColumn;
  // Holds the current column to which the next name entry will be read from.
  std::atomic<unsigned int> currentNameCol;
  // Holds the sample rate at which the data is automatically flushed to the file.
  unsigned int flushSampleRate{50};
  // Holds the minimum size for the data buffer, underload will trigger a file-read.
  unsigned int minBufferSize{20};
  // Holds the list of column names.
  std::vector<std::string> names;
  // Holds the map from the column names to the index within a value-row.
  mutable std::map<std::string, std::size_t> named_indices;
  // Holds the data buffer.
  std::queue<double> values_rm;
  // Holds the input-stream of the data record.
  std::shared_ptr<std::istream> in_stream;

  // Mutex to lock the read/write on the data buffer.
  std::mutex access_mutex;
  // Holds the instance of the data writing thread.
  std::shared_ptr<std::thread> reading_thread;

  /**
   * This class is used as a callable function-object for data writing thread.
   */
  struct extract_process {
   public:
    data_extractor& parent;
    explicit extract_process(data_extractor& aParent) : parent(aParent) {}
    void operator()();
  };

  void closeExtractProcess();

  /**
   * Overridable function which writes data to the file in whichever format specific to the derived class.
   */
  virtual bool readRow() { return true; }
  /**
   * Overridable function which writes column names to the file in whichever format specific to the derived class.
   */
  virtual bool readNames() { return true; }
  /**
   * Sets the filename for the file, overridable for the derived class which handles the file IO (or other).
   */
  virtual void setStreamImpl(
      const std::shared_ptr<std::istream>& aStreamPtr) = 0;

  void setStreamWrappedCall(const std::shared_ptr<std::istream>& aStreamPtr);

 public:
  /**
   * Returns the flushing sample rate of this data streamer in Hz.
   * \note A flushing sample rate of 0 signifies immediate transmission.
   * \return The flushing sample rate of this data streamer in Hz.
   */
  unsigned int getFlushSampleRate() const { return flushSampleRate; }

  /**
   * Sets the flushing sample rate of this data streamer in Hz.
   * \note A flushing sample rate of 0 signifies immediate transmission.
   * \param aFlushSampleRate The flushing sample rate of this data streamer in Hz.
   */
  void setFlushSampleRate(unsigned int aFlushSampleRate) {
    flushSampleRate = aFlushSampleRate;
  }

  /**
   * Returns the minimum size of the data buffer, after which, data must be replenished from the stream.
   * \return The minimum size of the data buffer, after which, data must be replenished from the stream.
   */
  unsigned int getMinBufferSize() const { return minBufferSize; }

  /**
   * Sets the minimum size of the data buffer, after which, data must be replenished from the stream.
   * \param aMinBufferSize The minimum size of the data buffer, after which, data must be replenished from the stream.
   */
  void setMinBufferSize(unsigned int aMinBufferSize) {
    minBufferSize = aMinBufferSize;
  }

  /**
   * This function is the factory to create named-value-row objects to represent a row of entries, addressable by name.
   * \return A fresh object that is ready to accept all the values of a row of entries to the data recorder.
   */
  named_value_row getFreshNamedValueRow() const {
    if (named_indices.size() != names.size()) {
      for (std::size_t i = 0; i < names.size(); ++i) {
        named_indices[names[i]] = i;
      }
    }
    return named_value_row(named_indices);
  }

  /**
   * Gets the current vector of column names.
   * \return The current vector of column names.
   */
  const std::vector<std::string>& getNames() const { return names; }

  /**
   * This function returns the number of columns in this extractor.
   * \return The number of columns in this extractor.
   */
  unsigned int getColCount() const { return colCount; }

  /// Data record-specific flags for special operations.
  enum flag {
    end_value_row,  ///< Ends the recording of the data entry for the current row.
    advance,        ///< Reads the data buffer from the file.
    close           ///< Closes the file and data buffer is flushed.
  };

  /**
   * Default Constructor.
   */
  data_extractor()
      : shared_object(), colCount(0), currentColumn(0), currentNameCol(0) {}

  /**
   * Destructor.
   */
  ~data_extractor() override;

  /**
   * Operator to extract a data entry.
   */
  data_extractor& operator>>(double& value);

  /**
   * Operator to extract a column name.
   */
  data_extractor& operator>>(std::string& name);

  /**
   * Operator to extract an entire row of named entries.
   */
  data_extractor& operator>>(named_value_row& values);

  /**
   * Operator to issue a special operation's flag.
   */
  data_extractor& operator>>(flag some_flag);

  /**
   * Operator to extract a vector of in-order (nameless) entries.
   */
  template <typename Vector>
  data_extractor& operator>>(Vector& values) {
    for (std::size_t i = 0; i < values.size(); ++i) {
      (*this) >> values[i];
    }
    return *this;
  }

  /**
   * Sets the stream.
   */
  void setStream(std::istream& aStream) {
    setStreamWrappedCall(
        std::shared_ptr<std::istream>(&aStream, null_deleter()));
  }

  /**
   * Sets the stream via a shared-pointer.
   */
  void setStream(const std::shared_ptr<std::istream>& aStreamPtr) {
    setStreamWrappedCall(aStreamPtr);
  }

  /**
   * Sets the filename for the file, overridable for the derived class which handles the file IO (or other).
   */
  virtual void setFileName(const std::string& aFileName);

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override;
  void load(serialization::iarchive& A, unsigned int /*Version*/) override;

  RK_RTTI_MAKE_ABSTRACT_1BASE(data_extractor, 0x81200001, 1, "data_extractor",
                              shared_object)
};

}  // namespace ReaK::recorder

#endif  // REAK_CORE_RECORDERS_DATA_RECORD_H_
