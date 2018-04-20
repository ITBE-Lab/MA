
#ifndef TEST_RANSAC_H
#define TEST_RANSAC_H

/* Function to find mean of the array elements.
 */
template<typename TP>
TP Mean(const TP arr[], int n )
{
	// Calculate sum of all elements.
	TP sum = 0;
	for( int i = 0; i < n; i++ )
		sum = sum + arr[i];
	return sum / n;
} // function

// Function to find mean absolute
// deviation of given elements.
template<typename TP>
TP meanAbsoluteDeviation(const TP arr[], int n)
{   
	// Calculate the sum of absolute
	// deviation about mean.
	TP absSum = 0;
	for (int i = 0; i < n; i++)
		absSum = absSum + abs(arr[i] - Mean(arr, n));

	// Return mean absolute deviation about mean.
	return absSum / n;
} // function

std::pair<double, double> run_ransac(
        const std::vector<double>& rvxValues,
        const std::vector<double>& rvyValues
    );

void export_ransac();

#endif