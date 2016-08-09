/* Basic example showing the usage of QuantLibAdjoint
*/

#include "run_example_2.hpp"
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/instruments/makevanillaswap.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>

#include <boost/make_shared.hpp>

#include <iostream>

using namespace QuantLib;
using std::vector;
using std::cout;
using std::endl;
using std::ios;

void example_2() {
	Date referenceDate(3, Aug, 2016);
	Settings::instance().evaluationDate() = referenceDate;
	Actual365Fixed dayCounter;

	// Example 1

	// These will be the X (independent) and Y (dependent) vectors
	vector<Rate> zeroRate;
	zeroRate.push_back(.02);
	zeroRate.push_back(.022);
	zeroRate.push_back(.026);
	zeroRate.push_back(.027);
	zeroRate.push_back(.029);
	vector<Date> datesZeroRate;
	datesZeroRate.push_back(Date(9, August, 2016));
	datesZeroRate.push_back(Date(9, August, 2017));
	datesZeroRate.push_back(Date(9, August, 2018));
	datesZeroRate.push_back(Date(9, August, 2019));
	datesZeroRate.push_back(Date(9, August, 2020));
	vector<Rate> swapNpv(1, 0.0);


	// Start taping with zeroRate as independent variable and set up flat zero curve
	cl::Independent(zeroRate);
	/*RelinkableHandle<YieldTermStructure> flatCurve(boost::make_shared<FlatForward>(referenceDate, zeroRate[0], dayCounter));
	flatCurve->enableExtrapolation();*/
	
	RelinkableHandle<YieldTermStructure> zeroCurve(boost::make_shared<ZeroCurve>(datesZeroRate, zeroRate, dayCounter, Calendar()));
	zeroCurve->enableExtrapolation();

	// Create and price swap
	Period swapTenor(5, Years);
	boost::shared_ptr<IborIndex> iborIndex = boost::make_shared<Euribor6M>(zeroCurve);
	Rate fixedRate = 0.03;
	Period forwardStart(0, Days);
	boost::shared_ptr<VanillaSwap> swap = MakeVanillaSwap(swapTenor, iborIndex, fixedRate, forwardStart).withNominal(100);
	swapNpv[0] = swap->NPV();

	// Stop taping and transfer operation sequence to function f (ultimately an AD function object)
	cl::tape_function<double> f(zeroRate, swapNpv);

	// Calculate d(swapNpv) / d(zero) with forward and reverse mode
	vector<double> dZ(1, 1.0);
	double forwardDeriv = f.Forward(1, dZ)[0];
	double reverseDeriv = f.Reverse(1, dZ)[0];

	// Calculate analytically the derivative
	Real derivative = 0.0;
	const Leg& fixedLeg = swap->fixedLeg();
	for (const auto& cf : fixedLeg) {
		Real amount = cf->amount();
		Time time = dayCounter.yearFraction(referenceDate, cf->date());
		DiscountFactor discount = zeroCurve->discount(time);
		derivative += amount * time * discount;
	}
	Time timeToStart = dayCounter.yearFraction(referenceDate, swap->startDate());
	Time timeToEnd = dayCounter.yearFraction(referenceDate, swap->maturityDate());
	derivative += 100 * (timeToEnd * zeroCurve->discount(timeToEnd) - timeToStart * zeroCurve->discount(timeToStart));

	//Compare bumped value to Taylor approximation
	Real basis_point = 0.01;
	Real approx_sensitivity = derivative * basis_point;
	//new bumped yield curve
	//when yield curve changes, swap recalculates 
	//flatCurve.linkTo(boost::make_shared<FlatForward>(referenceDate, zeroRate[0] + basis_point, dayCounter));
	//zeroCurve.linkTo(boost::make_shared<ZeroCurve>(datesZeroRate, zeroRate + basis_point, dayCounter));
	Real bumped_npv = swap->NPV();


	// Output the results
	cout.precision(9);
	cout.setf(ios::fixed, ios::floatfield);
	cout << "Zero Curve" << endl;
	cout << "Forward derivative:  " << forwardDeriv << endl;
	cout << "Reverse derivative:  " << reverseDeriv << endl;
	cout << "Analytic derivative: " << derivative << endl;



}