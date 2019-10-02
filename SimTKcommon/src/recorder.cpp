/*
 *      Recorder -- A package for Algorithmic Differentiation with CasADi
 *
 *      Copyright (C) 2019 The Authors
 *      Author: Joris Gillis
 *      Contributor: Antoine Falisse
 *
 *      Licensed under the Apache License, Version 2.0 (the "License"); you
 *      may not use this file except in compliance with the License. You may
 *      obtain a copy of the License at
 *      http://www.apache.org/licenses/LICENSE-2.0.
 *
 *      Unless required by applicable law or agreed to in writing, software
 *      distributed under the License is distributed on an "AS IS" BASIS,
 *      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 *      implied. See the License for the specific language governing
 *      permissions and limitations under the License.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include "SimTKcommon/internal/recorder.h"

static int counter_asserts = 0;

Recorder::Recorder(const Recorder& r) {
  value_ = r.value_;
  if (r.is_symbol()) {
    id_ = get_id();
    stream() << "a" + std::to_string(id_) << " = " << r.repr() << ";% copy constructor" << value_ << std::endl;
  } else {
    id_ = -1;
  }
}

Recorder::~Recorder() {
       if(id_ == 0) std::cout << "goodbye" << id_ << ":" << value_ << std::endl;
}

Recorder::Recorder() : id_(-1), value_(3.14) {}
Recorder::Recorder(double value) : id_(-1), value_(value) {}

void Recorder::operator<<=(double value) {
  if (is_symbol()) throw std::runtime_error("Needs to be symbolic");
  id_ = get_id();
  stream() << "if nom" << std::endl;
  stream() << "  " << repr() << " = " << value << ";" << std::endl;
  stream() << "else" << std::endl;
  stream() << "  " << repr() << " = x(" << counter_input+1 << ");" << std::endl;
  stream() << "end" << std::endl;
  counter_input++;
  value_ = value;
}

Recorder& Recorder::operator = ( const Recorder& r) {
  value_ = r.value_;
  if (r.is_symbol()) {
    id_ = get_id();
    stream() << "a" + std::to_string(id_) << " = " << r.repr() << ";% copy assignment" << value_ << std::endl;
  } else {
    id_ = -1;
  }
  return *this;
}

void Recorder::operator>>=(double& value) {
  if (!is_symbol()) throw std::runtime_error("Needs to be symbolic");
  stream() << "if ~nom" << std::endl;
  stream() << "y{" << counter_output+1 << "} = " << repr() << ";%" << value_  << std::endl;
  stream() << "end" << std::endl;
  counter_output++;
  value = value_;
}

double Recorder::getValue() const {return value_;}

Recorder operator+(const Recorder& lhs, const Recorder& rhs) {
  return Recorder::from_binary(lhs, rhs, lhs.value_ + rhs.value_, "plus");
}
Recorder operator*(const Recorder& lhs, const Recorder& rhs) {
  return Recorder::from_binary(lhs, rhs, lhs.value_ * rhs.value_, "times");
}
Recorder operator-(const Recorder& lhs, const Recorder& rhs) {
  return Recorder::from_binary(lhs, rhs, lhs.value_ - rhs.value_, "minus");
}
Recorder operator/(const Recorder& lhs, const Recorder& rhs) {
  return Recorder::from_binary(lhs, rhs, lhs.value_ / rhs.value_, "rdivide");
}
bool operator>=(const Recorder& lhs, const Recorder& rhs) {
  return static_cast<bool>(Recorder::from_binary(lhs, rhs, lhs.value_ >= rhs.value_, "ge"));
}
bool operator<=(const Recorder& lhs, const Recorder& rhs) {
  return static_cast<bool>(Recorder::from_binary(lhs, rhs, lhs.value_ <= rhs.value_, "le"));
}
bool operator>(const Recorder& lhs, const Recorder& rhs) {
  return static_cast<bool>(Recorder::from_binary(lhs, rhs, lhs.value_ > rhs.value_, "gt"));
}
bool operator<(const Recorder& lhs, const Recorder& rhs) {
  return static_cast<bool>(Recorder::from_binary(lhs, rhs, lhs.value_ < rhs.value_, "lt"));
}
bool operator!=(const Recorder& lhs, const Recorder& rhs) {
  return static_cast<bool>(Recorder::from_binary(lhs, rhs, lhs.value_ != rhs.value_, "ne"));
}
bool operator==(const Recorder& lhs, const Recorder& rhs) {
  return static_cast<bool>(Recorder::from_binary(lhs, rhs, lhs.value_ == rhs.value_, "eq"));
}
Recorder operator-(const Recorder& arg) {
    return Recorder::from_unary(arg, -arg.value_, "uminus");
}
Recorder pow( const Recorder&lhs, const Recorder& rhs) {
    return Recorder::from_binary(lhs, rhs, pow(lhs.value_,rhs.value_), "power");
}
Recorder fmax ( const Recorder&lhs, const Recorder& rhs) {
	return Recorder::from_binary(lhs, rhs, fmax(lhs.value_,rhs.value_), "max");
}
Recorder fmin ( const Recorder&lhs, const Recorder& rhs) {
	return Recorder::from_binary(lhs, rhs, fmin(lhs.value_,rhs.value_), "min");
}
Recorder atan2 ( const Recorder&lhs, const Recorder& rhs) {
	return Recorder::from_binary(lhs, rhs, atan2(lhs.value_,rhs.value_), "atan2");
}
Recorder exp(const Recorder& arg) {
    return Recorder::from_unary(arg, exp(arg.value_), "exp");
}
Recorder log(const Recorder& arg) {
    return Recorder::from_unary(arg, log(arg.value_), "log");
}
Recorder sqrt(const Recorder& arg) {
    return Recorder::from_unary(arg, sqrt(arg.value_), "sqrt");
}
Recorder sin(const Recorder& arg) {
    return Recorder::from_unary(arg, sin(arg.value_), "sin");
}
Recorder cos(const Recorder& arg) {
    return Recorder::from_unary(arg, cos(arg.value_), "cos");
}
Recorder tan(const Recorder& arg) {
    return Recorder::from_unary(arg, tan(arg.value_), "tan");
}
Recorder asin(const Recorder& arg) {
    return Recorder::from_unary(arg, asin(arg.value_), "asin");
}
Recorder acos(const Recorder& arg) {
    return Recorder::from_unary(arg, acos(arg.value_), "acos");
}
Recorder atan(const Recorder& arg) {
    return Recorder::from_unary(arg, atan(arg.value_), "atan");
}
Recorder log10(const Recorder& arg) {
    return Recorder::from_unary(arg, log10(arg.value_), "log10");
}
Recorder sinh(const Recorder& arg) {
    return Recorder::from_unary(arg, sinh(arg.value_), "sinh");
}
Recorder cosh(const Recorder& arg) {
    return Recorder::from_unary(arg, cosh(arg.value_), "cosh");
}
Recorder tanh(const Recorder& arg) {
    return Recorder::from_unary(arg, tanh(arg.value_), "tanh");
}
Recorder asinh(const Recorder& arg) {
    return Recorder::from_unary(arg, asinh(arg.value_), "asinh");
}
Recorder acosh(const Recorder& arg) {
    return Recorder::from_unary(arg, acosh(arg.value_), "acosh");
}
Recorder atanh(const Recorder& arg) {
    return Recorder::from_unary(arg, atanh(arg.value_), "atanh");
}
Recorder erf(const Recorder& arg) {
    return Recorder::from_unary(arg, erf(arg.value_), "erf");
}
Recorder fabs(const Recorder& arg) {
    return Recorder::from_unary(arg, fabs(arg.value_), "fabs");
}
Recorder ceil(const Recorder& arg) {
    return Recorder::from_unary(arg, ceil(arg.value_), "ceil");
}
Recorder floor(const Recorder& arg) {
    return Recorder::from_unary(arg, floor(arg.value_), "floor");
}

Recorder::operator bool() const {
   bool ret = value_==1;

    if (is_symbol()) {
      counter_bool++;
      stream() << "a{" << counter_asserts+1 << "} = " << repr() << "-" << value_ << ";%" << value_  << std::endl;
    }
    return ret;
}

std::ostream& operator<<(std::ostream &stream, const Recorder& obj) {
  obj.disp(stream);
  return stream;
}

std::istream& operator >> (std::istream& is, const Recorder& a) {
  throw std::runtime_error("No way!");
}

void Recorder::stop_recording() {
  stream() << "if ~nom, y = vertcat(y{:}); end;" << std::endl;
  if (counter_bool > 0) {
    stream() << "a = vertcat(a{:});" << std::endl;
  }
}

void Recorder::disp(std::ostream &stream) const {
  if (is_symbol()) {
    stream << "[#" << id_ << "|" << value_ << "]";
  } else {
    stream << "(" << value_ << ")";
  }
}

int Recorder::get_id() {
  counter++;
  return counter;
}

bool Recorder::is_symbol() const {
  return id_>=0;
}

std::string Recorder::repr() const {
  if (is_symbol()) {
    return "a" + std::to_string(id_);
  } else {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(16);
    ss << value_;
    return ss.str();
  }
}

bool is_suspicious(double v) {
  return (v>0 && v <1e-200) || (v<0 && v >-1e-200);
}

Recorder Recorder::from_binary(const Recorder& lhs, const Recorder& rhs, double res, const std::string& op) {
  if (lhs.is_symbol() || rhs.is_symbol()) {
    int id = get_id();
    stream() << "a" << id << " = " << op << "(" << lhs.repr() <<  "," <<   rhs.repr() << ");" << std::endl;
    stream() << "if nom, assert(" << "a" << id << "==" << res << "); end;" << std::endl;

    if (is_suspicious(res)) {
      stream() << "% suspicious activity" << std::endl;
    }

    return Recorder(res, id);
  } else {
    return Recorder(res);
  }
}

Recorder Recorder::from_unary(const Recorder& arg, double res, const std::string& op) {
  if (arg.is_symbol()) {
    int id = get_id();
    stream() << "a" << id << " = " << op << "(" << arg.repr() << ");" << std::endl;
    stream() << "if nom, assert(" << "a" << id << "==" << res << "); end;" << std::endl;

    if (is_suspicious(res)) {
      stream() << "% suspicious activity" << std::endl;
    }

    return Recorder(res, id);
  } else {
    return Recorder(res);
  }
}

class StreamWrapper {
public:
    StreamWrapper() {
        stream = new std::ofstream("foo.m");
        (*stream) << std::scientific << std::setprecision(16);
        (*stream) << "function [y,a]=foo(x)" << std::endl;
        (*stream) << "nom = nargin==0;" << std::endl;
    }

    std::ofstream* stream;
};
static StreamWrapper stream_wrapper{};

std::ofstream& Recorder::stream() {
  return *stream_wrapper.stream;
};

Recorder::Recorder(double value, int id) {
  id_ = id;
  value_ = value;
}

int Recorder::counter = 0;
int Recorder::counter_input = 0;
int Recorder::counter_output = 0;
int Recorder::counter_bool = 0;
