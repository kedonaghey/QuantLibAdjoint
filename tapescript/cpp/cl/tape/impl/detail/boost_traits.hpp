/*
Copyright (C) 2015-present CompatibL

Performance test results and finance-specific examples are available at:

http://www.modval.org/adjoint

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

// # include <boost/accumulators/framework/accumulator_base.hpp>

namespace boost
{
    namespace detail {
        template <typename Type> struct is_arithmetic_impl;

        template<typename Base>
        struct is_arithmetic_impl<cl::tape_wrapper<Base>>
        {
            typedef typename cl::remove_ad<typename cl::tape_wrapper<Base>::value_type>::type value_type;
#if defined BOOST_STATIC_CONSTANT
            BOOST_STATIC_CONSTANT(bool, value =
                (::boost::type_traits::ice_or<
                    ::boost::is_integral<value_type>::value,
                    ::boost::is_float<value_type>::value
                >::value));
#else
            static const bool value = std::is_arithmetic<value_type>::value;
#endif
            typedef typename std::is_arithmetic<value_type>::type type;
        };
    }
}

namespace cl
{
    template <class Base>
    class tape_wrapper;
}

namespace boost {  namespace accumulators
    {
        namespace impl
        {
            template<typename Sample>
            struct max_impl;

            ///////////////////////////////////////////////////////////////////////////////
            // max_impl
            template<typename Base>
            struct max_impl <cl::tape_wrapper<Base>>
                : accumulator_base
            {
                typedef cl::tape_wrapper<Base> Sample;
                // for boost::result_of
                typedef Sample result_type;

                template<typename Args>
                max_impl(Args const &args)
                    : max_(-std::numeric_limits<double>::max())
                {
#           if defined CL_COMPILE_TIME_DEBUG
#               pragma message ("boost overloading: " __FUNCSIG__)
#           endif
                }

                template<typename Args>
                void operator ()(Args const &args)
                {
#           if defined CL_COMPILE_TIME_DEBUG
#               pragma message ("boost overloading: " __FUNCSIG__)
#           endif

                    numeric::max_assign(this->max_, args[sample]);
                }

                template <class T>
                result_type result(T) const
                {
                    return this->max_;
                }

            private:
                Sample max_;
            };
        }
    }
}

