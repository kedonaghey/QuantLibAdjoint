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


#ifndef cl_adjoint_array_optimization_impl_hpp
#define cl_adjoint_array_optimization_impl_hpp
#pragma once


#include <ql/quantlib.hpp>
#include <cl/tape/impl/ql/blackcalculator_t.hpp>

#include "adjointarraytest.hpp"
#include "utilities.hpp"
#include "adjointtestutilities.hpp"
#include "adjointtestbase.hpp"

using namespace QuantLib;
using namespace boost::unit_test_framework;

#define BLACK_SCHOLES_FUNC // delta greek calculation test
//#define BLACK_SCHOLES_DELTA  // gamma greek calculation test
 
#define OUTPUT_FOLDER_NAME "AdjointArrayOptimization"

#if defined BLACK_SCHOLES_FUNC
#   define TARGET_FUNCTION(x) (black_func(x))
#   define TARGET_FUNCTION_NAME "BlackScholes price"
#elif defined BLACK_SCHOLES_DELTA
#   define TARGET_FUNCTION(x) (black_delta(x))
#   define TARGET_FUNCTION_NAME "BlackScholes delta"
#endif

namespace
{
    enum
    {
#if defined CL_GRAPH_GEN
        // Number of points for dependency graphics.
        pointNo = 50,
        // Number of points for performance graphic.
        iterNo = 50,
        // Step for portfolio size for performance testing .
        step = 10,
#else
        // Number of points for dependency graphics.
        pointNo = 0,
        // Number of points for performance graphic.
        iterNo = 0,
        // Step for portfolio size for performance testing .
        step = 1,
#endif
        // Defines performance accuracy. Its value is a minimum number
        // of calling of O(1) complexity methods per one performance test.
        iterNumFactor = 40000,
    };

    const double finite_diff_step = 1e-5;

    template <class T>
    inline T black_func(const T& x)
    {
        T spot = x;
        T strike = 100;
        T discount = 0.5;
        T forward = spot / discount;
        T stdDev = 0.2;

        return BlackCalculator_T<T>(Option::Call, strike, forward, stdDev, discount).value();
    }

    template <class T>
    inline T black_delta(const T& x)
    {
        T spot = x;
        T strike = 100;
        T discount = 0.5;
        T forward = spot / discount;
        T stdDev = 0.2;

        return BlackCalculator_T<T>(Option::Call, strike, forward, stdDev, discount).delta(x);
    }

    // Struct for graphics recording.
    struct OutStruct
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "X", "Adjoint/Tape results", "Fin. diff./Function results"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type& operator<<(stream_type& stm, OutStruct& v)
        {
            stm << v.x_
                << ";" << v.adj_
                << ";" << v.thr_
                << std::endl;
            return stm;
        }

        Real x_;
        Real adj_;
        Real thr_;
    };

    struct NoOptTest
        : cl::AdjointTest<NoOptTest>
    {
        static const CalcMethod default_method = reverse;

        explicit NoOptTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(size)
            , output_(size)
            , X_(size)
            , Y_()
            , analyticalY_(size)
        {
            setLogger(logger);
            for (size_t i = 0; i < size_; i++)
            {
                input_[i] = 50 + 10 * (double(i) / size_);
            }

            std::copy(input_.begin(), input_.end(), X_.begin());
        }

        void recordTape()
        {
            cl::Independent(X_);
            calcY();
            f_ = std::make_unique<cl::tape_function<double>>(X_, Y_);
        }

        void calcY()
        {
            std::transform(X_.begin(), X_.end(), output_.begin(), targetFunc);

            Real sum = std::accumulate(output_.begin(), output_.end(), Real());

            Y_ = { sum };
        }

        inline static Real targetFunc(const Real& x)
        {
            return TARGET_FUNCTION(x);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            Real h = finite_diff_step;
            analyticalResults_.resize(size_);
            analyticalY_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                Real left = targetFunc(input_[i] - h);
                Real right = targetFunc(input_[i] + h);
                analyticalResults_[i] = (right - left) / (2 * h);
                analyticalY_[i] = (right + left) / 2;
            }
        }

        Size indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return 2 * iterNumFactor; }

        double relativeTol() const { return 0.05; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<double> input_;
        std::vector<Real> output_;
        std::vector<Real> X_;
        std::vector<Real> Y_;
        std::vector<Real> analyticalY_;
    };
    
    struct NoOptTestData
    {
        explicit NoOptTestData()
        : outPerform_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " differentiation performance without optimization." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME, {
                { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance without optimization." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME, {
                { "title", "Tape size dependence on number of variables" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , outDeriv_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME "//output", {
                { "filename", "Derivative" }
            , { "not_clear", "Not" }
            , { "title", "Function derivative" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
            , outValue_(OUTPUT_FOLDER_NAME "//NoOpt " TARGET_FUNCTION_NAME "//output", {
                { "filename", "FunctionValue" }
            , { "not_clear", "Not" }
            , { "title", "Function value" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<NoOptTest> getTest(size_t size)
        {
            return std::make_shared<NoOptTest>(size, &outPerform_);
        }

        void recordDependencePlot(size_t pointNo)
        {
            std::vector<OutStruct> outDeriv(pointNo);
            std::vector<OutStruct> outValue(pointNo);
            auto test = getTest(pointNo);
            test->calcWithReverse();
            test->calcAnalytical();
            for (Size i = 0; i < pointNo; i++)
            {
                outDeriv[i] = { test->input_[i], test->reverseResults_[i], test->analyticalResults_[i] };
                outValue[i] = { test->input_[i], test->output_[i], test->analyticalY_[i] };
            }
            outDeriv_ << outDeriv;
            outValue_ << outValue;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outDeriv_;
        cl::tape_empty_test_output outValue_;
    };
    
    struct InnerArrayTest
        : cl::AdjointTest<InnerArrayTest>
    {
        typedef cl::tape_value inner_type;
        typedef cl::tape_wrapper<inner_type> tape_type;
        typedef cl::tape_function<inner_type> func_type;
        static const CalcMethod default_method = reverse;

        explicit InnerArrayTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(size)
            , output_(size)
            , X_()
            , Y_()
            , analyticalY_(size)
            , adf_()
        {
            setLogger(logger);
            for (size_t i = 0; i < size_; i++)
            {
                input_[i] = 50 + 10 * (double(i) / size_);
            }

            tape_type x(inner_type(std::valarray<double>(input_.data(), size)));
            X_ = { std::move(x) };
        }

        void recordTape()
        {
            cl::Independent(X_);
            calcY();
            adf_ = std::make_unique<func_type>(X_, Y_);
            auto& res = CppAD::Value(Y_[0].value()).array_value_;
            output_.assign(std::begin(res), std::end(res));
        }

        void calcY()
        {
            Y_ = { targetFunc(X_[0]) };
        }

        template <class T>
        inline static T targetFunc(const T& x)
        {
            return TARGET_FUNCTION(x);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = finite_diff_step;
            analyticalResults_.resize(size_);
            analyticalY_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                double left = targetFunc(input_[i] - h);
                double right = targetFunc(input_[i] + h);
                analyticalResults_[i] = (right - left) / (2 * h);
                analyticalY_[i] = (right + left) / 2;
            }
        }

        void calcForward()
        {
            typedef std::vector<inner_type> Vec;
            Vec dy = adf_->Forward(1, Vec{ 1 });
            setForwardResults(dy[0]);
        }

        void calcReverse()
        {
            typedef std::vector<inner_type> Vec;
            Vec dy = adf_->Reverse(1, Vec{ 1 });
            setReverseResults(dy[0]);
        }

        void setForwardResults(const inner_type& res)
        {
            auto& result = res.array_value_;
            forwardResults_.assign(std::begin(result), std::end(result));            
        }

        void setReverseResults(const inner_type& res)
        {
            auto& result = res.array_value_;
            reverseResults_.assign(std::begin(result), std::end(result));
        }

        size_t memory()
        {
            if (!adf_)
            {
                return 0;
            }
            return adf_->Memory();
        }

        Size indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return 5 * iterNumFactor; }

        size_t tapePerfIteration()
        {
            return iterNumFactor / size_ + 1;
        }

        double relativeTol() const { return 0.05; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<double> input_;
        std::vector<double> output_;
        std::vector<tape_type> X_;
        std::vector<tape_type> Y_;
        std::vector<double> analyticalY_;
        std::unique_ptr<func_type> adf_;
    };
    
    struct InnerArrayTestData
    {
        explicit InnerArrayTestData()
        : outPerform_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " differentiation performance using InnerArray." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME, {
                { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance using InnerArray." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME, {
                { "title", "Tape size dependence on number of variables"
                "\\n (memory dynamically allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , outDeriv_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME "//output", {
                { "filename", "Derivative" }
            , { "not_clear", "Not" }
            , { "title", "Function derivative" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
            , outValue_(OUTPUT_FOLDER_NAME "//InnerArray " TARGET_FUNCTION_NAME "//output", {
                { "filename", "FunctionValue" }
            , { "not_clear", "Not" }
            , { "title", "Function value" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            outPerform_.log() << "Start InnerArrayTestData output.\n";
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<InnerArrayTest> getTest(size_t size)
        {
            return std::make_shared<InnerArrayTest>(size, &outPerform_);
        }

        void recordDependencePlot(size_t pointNo)
        {
            std::vector<OutStruct> outDeriv(pointNo);
            std::vector<OutStruct> outValue(pointNo);
            auto test = getTest(pointNo);
            test->recordTape();
            test->calcReverse();
            test->calcAnalytical();
            for (Size i = 0; i < pointNo; i++)
            {
                outDeriv[i] = { test->input_[i], test->reverseResults_[i], test->analyticalResults_[i] };
                outValue[i] = { test->input_[i], test->output_[i], test->analyticalY_[i] };
            }
            outDeriv_ << outDeriv;
            outValue_ << outValue;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outDeriv_;
        cl::tape_empty_test_output outValue_;
    };
    
    struct MixedTest
        : cl::AdjointTest<MixedTest>
    {
        typedef cl::tape_value inner_type;
        typedef cl::tape_wrapper<inner_type> tape_type;
        typedef cl::tape_function<inner_type> func_type;
        static const CalcMethod default_method = reverse;

        explicit MixedTest(Size size, cl::tape_empty_test_output* logger = nullptr)
            : AdjointTest()
            , size_(size)
            , input_(size)
            , output_(size)
            , X_()
            , Y_()
            , analyticalY_(size)
            , adf_()
            , x_val_()
        {
            setLogger(logger);
            for (size_t i = 0; i < size_; i++)
            {
                input_[i] = 50 + 10 * (double(i) / size_);
            }

            x_val_ = std::valarray<double>(input_.data(), size);
            X_ = { tape_type(50) };
        }

        void recordTape()
        {
            cl::Independent(X_);
            calcY();
            adf_ = std::make_unique<func_type>();
            adf_->Dependent(X_, Y_);
            typedef std::vector<inner_type> Vec;
            auto res = adf_->Forward(0, Vec{ x_val_ })[0].array_value_;
            output_.assign(std::begin(res), std::end(res));
        }

        void calcY()
        {
            Y_ = { targetFunc(X_[0]) };
        }

        template <class T>
        inline static T targetFunc(const T& x)
        {
            return TARGET_FUNCTION(x);
        }

        // Calculates derivatives using finite difference method.
        void calcAnalytical()
        {
            double h = finite_diff_step;
            analyticalResults_.resize(size_);
            analyticalY_.resize(size_);
            for (Size i = 0; i < size_; i++)
            {
                double left = targetFunc(input_[i] - h);
                double right = targetFunc(input_[i] + h);
                analyticalResults_[i] = (right - left) / (2 * h);
                analyticalY_[i] = (right + left) / 2;
            }
        }

        void calcForward()
        {
            typedef std::vector<inner_type> Vec;
            Vec dy = adf_->Forward(1, Vec{ 1 });
            setForwardResults(dy[0]);
        }

        void calcReverse()
        {
            typedef std::vector<inner_type> Vec;
            Vec dy = adf_->Reverse(1, Vec{ 1 });
            setReverseResults(dy[0]);
        }

        void setForwardResults(const inner_type& res)
        {
            auto& result = res.array_value_;
            forwardResults_.assign(std::begin(result), std::end(result));
        }

        void setReverseResults(const inner_type& res)
        {
            auto& result = res.array_value_;
            reverseResults_.assign(std::begin(result), std::end(result));
        }

        size_t memory()
        {
            if (!adf_)
            {
                return 0;
            }
            return adf_->Memory();
        }

        Size indepVarNumber()
        {
            return size_;
        }

        size_t minPerfIteration() { return 10 * iterNumFactor; }

        size_t tapePerfIteration()
        {
            return iterNumFactor / size_ + 1;
        }

        double relativeTol() const { return 0.05; }

        double absTol() const { return 1e-5; }

        Size size_;
        std::vector<double> input_;
        std::vector<double> output_;
        std::vector<tape_type> X_;
        std::vector<tape_type> Y_;
        std::vector<double> analyticalY_;
        std::unique_ptr<func_type> adf_;
        inner_type x_val_;
    };
    
    struct MixedTestData
    {
        explicit MixedTestData()
        : outPerform_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " differentiation performance with Mixed optimization." }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "line_box_width", "-5" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        })
            , outAdjoint_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME, {
                { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance with Mixed optimization." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        })
            , outSize_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME, {
                { "title", "Tape size dependence on number of variables"
                "\\n (memory dynamically allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (MB)" }
            , { "cleanlog", "false" }
        })
            , outDeriv_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME "//output", {
                { "filename", "Derivative" }
            , { "not_clear", "Not" }
            , { "title", "Function derivative" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
            , outValue_(OUTPUT_FOLDER_NAME "//Mixed " TARGET_FUNCTION_NAME "//output", {
                { "filename", "FunctionValue" }
            , { "not_clear", "Not" }
            , { "title", "Function value" }
            , { "ylabel", "Y" }
            , { "xlabel", "X" }
            , { "cleanlog", "false" }
        })
        { }

        bool makeOutput()
        {
            outPerform_.log() << "Start MixedTestData output.\n";
            bool ok = true;
            if (pointNo > 0)
            {
                recordDependencePlot(pointNo);
            }
            ok &= cl::recordPerformance(*this, iterNo, step);
            return ok;
        }

        std::shared_ptr<MixedTest> getTest(size_t size)
        {
            return std::make_shared<MixedTest>(size, &outPerform_);
        }

        void recordDependencePlot(size_t pointNo)
        {
            std::vector<OutStruct> outDeriv(pointNo);
            std::vector<OutStruct> outValue(pointNo);
            auto test = getTest(pointNo);
            test->recordTape();
            test->calcReverse();
            test->calcAnalytical();
            for (Size i = 0; i < pointNo; i++)
            {
                outDeriv[i] = { test->input_[i], test->reverseResults_[i], test->analyticalResults_[i] };
                outValue[i] = { test->input_[i], test->output_[i], test->analyticalY_[i] };
            }
            outDeriv_ << outDeriv;
            outValue_ << outValue;
        }

        cl::tape_empty_test_output outPerform_;
        cl::tape_empty_test_output outAdjoint_;
        cl::tape_empty_test_output outSize_;
        cl::tape_empty_test_output outDeriv_;
        cl::tape_empty_test_output outValue_;
    };

    struct MethodsCompare
    {
        static std::deque<std::string > get_columns()
        {
            static std::deque<std::string > columns =
            {
                "N", "No Optimization", "InnerArray", "Mixed"
                , "Finite diff. (double)", "Finite diff. (Real)", "Finite diff. (InnerArray)"
            };

            return columns;
        }

        template <typename stream_type>
        friend inline stream_type& operator<<(stream_type& stm, MethodsCompare& v)
        {
            stm << v.n_
                << ";" << v.NoOpt_
                << ";" << v.InnerArrayOpt_
                << ";" << v.Mixed_
                << ";" << v.double_
                << ";" << v.Real_
                << ";" << v.FinInner_
                << std::endl;
            return stm;
        }

        Real n_;
        Real NoOpt_;
        Real InnerArrayOpt_;
        Real Mixed_;
        Real double_;
        Real Real_;
        Real FinInner_;
    };

    void compareSize()
    {
        cl::tape_empty_test_output outTotal(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation and tape recording"
            "\\n and function construction (including Forward(0) sweep) performance." }
            , { "filename", "Total" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "true" }
        });
        cl::tape_empty_test_output outTotalSmooth(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation and tape recording"
            "\\n and function construction (including Forward(0) sweep) performance." }
            , { "filename", "TotalSmooth" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "true" }
            , { "smooth", "default" }
        });
        cl::tape_empty_test_output outAdjoint(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance." }
            , { "filename", "Adjoint" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
        });
        cl::tape_empty_test_output outAdjointSmooth(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " Adjoint differentiation performance." }
            , { "filename", "AdjointSmooth" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            , { "smooth", "default" }
        });
        cl::tape_empty_test_output outTape(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", TARGET_FUNCTION_NAME " tape recording performance." }
            , { "filename", "Tape" }
            , { "not_clear", "Not" }
            , { "ylabel", "Time (s)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
            //, { "smooth", "default" }
        });
        cl::tape_empty_test_output outSize(OUTPUT_FOLDER_NAME "//Compare " TARGET_FUNCTION_NAME, {
            { "title", "Tape size dependence on number of variables"
            "\\n (memory dynamically allocated for arrays is not included)" }
            , { "filename", "TapeSize" }
            , { "not_clear", "Not" }
            , { "ylabel", "Memory (KB)" }
            , { "xlabel", "Number of variables" }
            , { "cleanlog", "false" }
        });

        size_t graphSize = iterNo + step - 1;
        std::vector<MethodsCompare> total(graphSize)
            , adjoint(graphSize)
            , tape(graphSize)
            , mem(graphSize);
        for (size_t i = 0; i < graphSize; i++)
        {
            size_t n;
            if (i < step)
            {
                n = i + 1;
            }
            else
            {
                n = (i - step + 2) * step;
            }
            auto no = NoOptTest(n, &outTotal);
            auto inn = InnerArrayTest(n, &outTotal);
            auto mix = MixedTest(n, &outTotal);
            auto chp = MixedTest(n, &outTotal);

            tape[i] = {
                Real(n)
                , no.testTapePerformance()
                , inn.testTapePerformance()
                , mix.testTapePerformance()
                , 0
                , 0
                , 0
            };
            adjoint[i] = {
                Real(n)
                , no.testPerformance<>()
                , inn.testPerformance<>()
                , mix.testPerformance<>()
                , inn.testPerformance<InnerArrayTest::analytical>()
                , no.testPerformance<NoOptTest::analytical>()
                , mix.testPerformance<MixedTest::analytical>()
            };
            total[i] = {
                Real(n)
                , tape[i].NoOpt_            + adjoint[i].NoOpt_
                , tape[i].InnerArrayOpt_    + adjoint[i].InnerArrayOpt_
                , tape[i].Mixed_            + adjoint[i].Mixed_
                , tape[i].double_           + adjoint[i].double_
                , tape[i].Real_             + adjoint[i].Real_
                , tape[i].FinInner_         + adjoint[i].FinInner_
            };
            mem[i] = {
                Real(n)
                , no.memory() / 1024.0
                , inn.memory() / 1024.0
                , mix.memory() / 1024.0
                , 0
                , 0
                , 0
            };
        }

        outTotal << total;
        outAdjoint << adjoint;
        outTotalSmooth << total;
        outAdjointSmooth << adjoint;
        outTape << tape;
        outSize << mem;
    }
}

#endif // #ifndef cl_adjoint_array_optimization_impl_hpp
