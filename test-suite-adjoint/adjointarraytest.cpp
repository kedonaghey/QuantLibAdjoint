/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 CompatibL

This file is part of QuantLib, a free-software/open-source library
for financial quantitative analysts and developers - http://quantlib.org/

QuantLib is free software: you can redistribute it and/or modify it
under the terms of the QuantLib license.  You should have received a
copy of the license along with this program; if not, please email
<quantlib-dev@lists.sf.net>. The license is also available online at
<http://quantlib.org/license.shtml>.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#define CL_BASE_SERIALIZER_OPEN
#define CL_TAPE_INNER_ARRAY_ENABLED

#include "adjointarrayimpl.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;


bool AdjointArrayTest::testNoOpt()
{
    BOOST_TEST_MESSAGE("Testing Adjoint for vector without optimization...");

    // Initializations...
    NoOptTestData td;
    Size n = 5;
    NoOptTest test(n, &td.outPerform_);
    std::vector<Real>& X = test.X_;

    // Declare input vector as independent variables.
    cl::Independent(X);

    // Calculate the sum of all function values.
    Real sum = 0.0;
    for (Size i = 0; i < n; i++)
    {
        sum += TARGET_FUNCTION(X[i]);
    }
    // Use the sum as one output variable.
    std::vector<Real> Y = { sum };

    // Declare tape function and stop the tape recording.
    cl::tape_function<double> f(X, Y);

    // Reverse mode calculations.
    test.reverseResults_ = f.Reverse(1, std::vector<double>{ 1 });

    // Result checking and output generation...
    test.calcAnalytical();

    std::cout << "Testing Adjoint for vector without optimization..." << std::endl;
    std::cout << "Adjoint derivatives\tFinite-difference derivatives" << std::endl;
    for (int i = 0; i < test.size_; i++)
    {
        std::cout << std::setprecision(8) << test.reverseResults_[i] << "\t\t" << test.analyticalResults_[i] << std::endl;
    }

    bool ok = test.checkAdjoint<>();
    ok &= test.testAdjoint<>();
    ok &= td.makeOutput();
    return ok;
}

bool AdjointArrayTest::testInnerArray()
{
    BOOST_TEST_MESSAGE("Testing Adjoint using InnerArray class...");

    // Initializations...
    InnerArrayTestData td;
    Size n = 5;
    InnerArrayTest test(n, &td.outPerform_);
    std::vector<cl::tape_object> X = { test.X_[0] };

    // Declare an independent variable.
    cl::Independent(X);

    // Use only one function call to calculate the whole vector.
    std::vector<cl::tape_object> Y = { TARGET_FUNCTION(X[0]) };

    // Declare a tape function and stop the tape recording.
    cl::tape_function<cl::tape_value> f(X, Y);

    // Forward mode calculation.
    cl::tape_value forw = f.Forward(1, std::vector<cl::tape_value>{ 1.0 })[0];

    // Reverse mode calculation.
    cl::tape_value rev = f.Reverse(1, std::vector<cl::tape_value>{ 1.0 })[0];

    // Result checking and output generation...
    test.setForwardResults(forw);
    test.setReverseResults(rev);
    test.calcAnalytical();

    std::cout << "Testing Adjoint using InnerArray class..." << std::endl;
    std::cout << "Adjoint derivatives\tFinite-difference derivatives" << std::endl;
    for (int i = 0; i < test.size_; i++)
    {
        std::cout << std::setprecision(8) << test.reverseResults_[i] << "\t\t" << test.analyticalResults_[i] << std::endl;
    }

    bool ok = test.check();
    ok &= test.test();
    ok &= td.makeOutput();
    return ok;
}

bool AdjointArrayTest::testMixed()
{
    BOOST_TEST_MESSAGE("Testing Adjoint using mixed optimization...");
    
    // Initializations...
    MixedTestData td;
    Size n = 5;
    MixedTest test(n, &td.outPerform_);
    cl::tape_value& input = test.x_val_;

    // Set the scalar value (default) as input variable.
    std::vector<cl::tape_object> X = { cl::tape_object(X_DEFAULT) };

    // Declare an independent variable.
    cl::Independent(X);

    // Calculate the target function for scalar valued input (during tape recording).
    std::vector<cl::tape_object> Y = { TARGET_FUNCTION(X[0]) };

    // Use default constructor for tape function (tape is still recorded).
    cl::tape_function<cl::tape_value> f;
    // Declare the dependent variable and stop the tape recording.
    f.Dependent(X, Y);
    // Initial Forward(0) sweep for actual data.
    f.Forward(0, std::vector<cl::tape_value>{ input });    

    // Forward mode calculation.
    cl::tape_value forw = f.Forward(1, std::vector<cl::tape_value>{ 1.0 })[0];

    // Reverse mode calculation.
    cl::tape_value rev = f.Reverse(1, std::vector<cl::tape_value>{ 1.0 })[0];

    // Result checking and output generation...
    test.setForwardResults(forw);
    test.setReverseResults(rev);
    test.calcAnalytical();

    std::cout << "Testing Adjoint using mixed optimization..." << std::endl;
    std::cout << "Adjoint derivatives\tFinite-difference derivatives" << std::endl;
    for (int i = 0; i < test.size_; i++)
    {
        std::cout << std::setprecision(8) << test.reverseResults_[i] << "\t\t" << test.analyticalResults_[i] << std::endl;
    }

    bool ok = test.check();
    ok &= test.test();
    ok &= td.makeOutput();
    return ok;
}

bool AdjointArrayTest::printAll()
{
    BOOST_TEST_MESSAGE("Testing Adjoint (comparing different methods)...");

#if defined CL_GRAPH_GEN
    compareSize();
#endif

    return true;
}

test_suite* AdjointArrayTest::suite()
{
    test_suite* suite = BOOST_TEST_SUITE("Adjoint array tests.");

    suite->add(QUANTLIB_TEST_CASE(&AdjointArrayTest::testNoOpt));
    suite->add(QUANTLIB_TEST_CASE(&AdjointArrayTest::testInnerArray));
    suite->add(QUANTLIB_TEST_CASE(&AdjointArrayTest::testMixed));
    suite->add(QUANTLIB_TEST_CASE(&AdjointArrayTest::printAll));

    return suite;
}


#ifdef CL_ENABLE_BOOST_TEST_ADAPTER

BOOST_AUTO_TEST_SUITE(ad_array)

BOOST_AUTO_TEST_CASE(testArrayNoOpt)
{
    BOOST_CHECK(AdjointArrayTest::testNoOpt());
}
BOOST_AUTO_TEST_CASE(testArrayInnerArray)
{
    BOOST_CHECK(AdjointArrayTest::testInnerArray());
}
BOOST_AUTO_TEST_CASE(testArrayMixed)
{
    BOOST_CHECK(AdjointArrayTest::testMixed());
}
BOOST_AUTO_TEST_CASE(testArrayAll)
{
    BOOST_CHECK(AdjointArrayTest::printAll());
}

BOOST_AUTO_TEST_SUITE_END()

#endif
