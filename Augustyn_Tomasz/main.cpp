#include "opencv2/opencv.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <windows.h>


using namespace cv;
using namespace std;

#define M_PI           3.14159265358979323846  /* pi */

//declarations
struct AngleAndElemNumber;
struct DistanceAndElemNumber;
struct AreaAndElemNumber;
enum JellyColors;
Point2f CalculateCenterPointOfPoints( std::vector<Point2f> mc2 );
std::vector<Point2f> RelateVectorOfPointsToTheCenterPoint( std::vector<Point2f> mc2, Point2f oCenterPoint );
std::vector<AngleAndElemNumber> ConvertToPolarCoordinates( std::vector<Point2f> mc2_related );
std::vector<DistanceAndElemNumber> CalculateDistanceToCenter( std::vector<Point2f> mc2, cv::Mat frame );
int Partition( std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r );
int Partition( std::vector<DistanceAndElemNumber> &VectorOfDistAndElem, int p, int r );
int Partition(std::vector<DistanceAndElemNumber> &VectorOfAreaAndElem, int p, int r);
void Quicksort( std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r );
void Quicksort( std::vector<DistanceAndElemNumber> &VectorOfDistAndElem, int p, int r );
void Quicksort(std::vector<DistanceAndElemNumber> &VectorOfAreaAndElem, int p, int r);
std::vector<Point2f> RemoveContoursFarestFromCenter( std::vector<DistanceAndElemNumber> VectorOfDistAndElem, std::vector<Point2f> mc2, std::vector< vector <Point> > contours_reduced_2, std::vector< vector <Point> > &contours_reduced_3 );
std::vector<Point2f> EnumerateVerticiesClockwise( std::vector<AngleAndElemNumber> VectorOfAngleAndElem, std::vector<Point2f> mc2 );
std::vector< vector <Point> > RemoveSmallAndBigContours( std::vector< vector <Point> > contours, double dLowerAreaThreshold, double dHigherAreaThreshold );
std::vector< vector <Point> > RemoveObjEnclosingCircleLessThen( std::vector< vector <Point> > contours_reduced, double MinimumRadius );
std::vector<Point2f> GetContoursMassCenters(std::vector< vector <Point> > contours_reduced_2);
std::vector< vector <Point> > ApproximateContours( std::vector< vector <Point> > contours_reduced_2, double precision );
double GetAverageNrOfPointsForContoursSet( std::vector< vector <Point> > contours_approxed );
std::vector<Point2f> GetTargetPointsDependingOnMarkerType( std::vector< vector <Point> > contours_approxed );
std::vector<AreaAndElemNumber> CalculateContoursArea( std::vector< vector <Point> > contours_jelly_reduced );
double Median( std::vector<AreaAndElemNumber> VectorOfAreaAndElem );
int Median( std::vector<int> color_jellies );
int CountJellies( std::vector<AreaAndElemNumber> VectorOfAreaAndElem );
int CountColorJelly( cv::Mat warp_hsv, int color );

//structure used in sorting algorithm (quicksort)
struct AngleAndElemNumber
{
	double dPolarAngle;
	int iElementNumber;
};

//structure used in sorting algorithm (quicksort)
struct DistanceAndElemNumber
{
	double dDistanceToCenter;
	int iElementNumber;
};

//structure used in sorting algorithm (quicksort)
struct AreaAndElemNumber
{
	double dArea;
	int iElementNumber;
};

enum JellyColors
{
	WHITE = 0,
	YELLOW = 1,
	ORANGE = 2,
	LIGHT_RED = 3,
	DARK_RED = 4,
	GREEN = 5
};

Point2f CalculateCenterPointOfPoints( std::vector<Point2f> mc3 )
{
	size_t iNumberOfPoints = mc3.size();
	float dCenterX = 0, dCenterY = 0;
	Point2f oCenterPoint;

	for (size_t i = 0; i < iNumberOfPoints; i++)
	{
		dCenterX = dCenterX + mc3[i].x;
		dCenterY = dCenterY + mc3[i].y;
	}
	dCenterX = dCenterX / iNumberOfPoints;
	dCenterY = dCenterY / iNumberOfPoints;
	oCenterPoint.x = dCenterX;
	oCenterPoint.y = dCenterY;

	return oCenterPoint;

}

std::vector<Point2f> RelateVectorOfPointsToTheCenterPoint( std::vector<Point2f> mc3, Point2f oCenterPoint )
{
	std::vector<Point2f> mc3_related;
	size_t iNumberOfPoints = mc3.size();
	Point2f temp_point;

	for (size_t i = 0; i < iNumberOfPoints; i++)
	{
		temp_point.x = mc3[i].x - oCenterPoint.x;
		temp_point.y = mc3[i].y - oCenterPoint.y;
		mc3_related.push_back(temp_point);
	}
	return mc3_related;

}

std::vector<AngleAndElemNumber> ConvertToPolarCoordinates( std::vector<Point2f> mc3_related )
{
	size_t iNumberOfPoints = mc3_related.size();
	std::vector<AngleAndElemNumber> VectorOfAngleAndElem;
	AngleAndElemNumber TempStruct;
	float tempX = 0.0, tempY = 0.0;

	for ( size_t i = 0; i < iNumberOfPoints; i++ )
	{
		TempStruct.iElementNumber = (int)i;
		tempX = mc3_related[i].x;
		tempY = mc3_related[i].y;

		if ( ( tempX > 0 ) && ( tempY >= 0 ) )
		{
			TempStruct.dPolarAngle = atan(tempY / tempX);
		}
		else if ( ( tempX > 0) && (tempY < 0 ) )
		{
			TempStruct.dPolarAngle = atan(tempY / tempX) + 2 * M_PI;
		}
		else if ( tempX < 0 )
		{
			TempStruct.dPolarAngle = atan(tempY / tempX) + M_PI;
		}
		else if ( (tempX == 0) && (tempY > 0) )
		{
			TempStruct.dPolarAngle = M_PI / 2;
		}
		else if ( (tempX == 0) && (tempY < 0) )
		{
			TempStruct.dPolarAngle = 3 * M_PI / 2;
		}
		else
		{
			TempStruct.dPolarAngle = 0;
		}

		VectorOfAngleAndElem.push_back( TempStruct );
	}

	return VectorOfAngleAndElem;
}

std::vector<DistanceAndElemNumber> CalculateDistanceToCenter( std::vector<Point2f> mc2, cv::Mat frame )
{
	size_t iNumberOfPoints = mc2.size();
	std::vector<DistanceAndElemNumber> VectorOfDistAndElem;
	DistanceAndElemNumber TempStruct;
	double tempX = 0.0, tempY = 0.0;
	double CenterX = frame.cols / 2;
	double CenterY = frame.rows / 2;
	double distance = 0.0;
	

	for (size_t i = 0; i < iNumberOfPoints; i++)
	{
		TempStruct.iElementNumber = (int)i;
		tempX = mc2[i].x;
		tempY = mc2[i].y;
		TempStruct.dDistanceToCenter = sqrt(pow(CenterX - tempX, 2) + pow(CenterY - tempY, 2));

		VectorOfDistAndElem.push_back(TempStruct);
	}

	return VectorOfDistAndElem;
}


//Overloaded partition function
int Partition( std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r ) //we divide structure in 2 parts, in the first one all angle degrees are less or equal to x, in the second one greater or equal to x
{
	double x = VectorOfAngleAndElem[p].dPolarAngle; // we take x
	int i = p, j = r, index; // i, j - indexes in structure
	double w;
	while (true) // infinite loop - we leave it only by returning j
	{
		while (VectorOfAngleAndElem[j].dPolarAngle > x) // when elements are greater than x 
			j--;
		while (VectorOfAngleAndElem[i].dPolarAngle < x) // when elements are less than x 
			i++;
		if (i < j) // we swap places when i < j
		{
			w = VectorOfAngleAndElem[i].dPolarAngle;
			index = VectorOfAngleAndElem[i].iElementNumber;

			VectorOfAngleAndElem[i].dPolarAngle = VectorOfAngleAndElem[j].dPolarAngle;
			VectorOfAngleAndElem[i].iElementNumber = VectorOfAngleAndElem[j].iElementNumber;

			VectorOfAngleAndElem[j].dPolarAngle = w;
			VectorOfAngleAndElem[j].iElementNumber = index;
			i++;
			j--;
		}
		else // when i >= j  we return j as the structure division point
			return j;
	}
}

//Overloaded partition function
int Partition(std::vector<DistanceAndElemNumber> &VectorOfDistAndElem, int p, int r) // we divide structure in 2 parts, in the first one all distances are less or equal to x, in the second one greater or equal to x
{
	double x = VectorOfDistAndElem[p].dDistanceToCenter; // we take x
	int i = p, j = r, index; // i, j - indexes in structure
	double w;
	while (true) // infinite loop - we leave it only by returning j
	{
		while (VectorOfDistAndElem[j].dDistanceToCenter > x) // when elements are greater than x 
			j--;
		while (VectorOfDistAndElem[i].dDistanceToCenter < x) // when elements are less than x
			i++;
		if (i < j) // we swap places when i < j
		{
			w = VectorOfDistAndElem[i].dDistanceToCenter;
			index = VectorOfDistAndElem[i].iElementNumber;

			VectorOfDistAndElem[i].dDistanceToCenter = VectorOfDistAndElem[j].dDistanceToCenter;
			VectorOfDistAndElem[i].iElementNumber = VectorOfDistAndElem[j].iElementNumber;

			VectorOfDistAndElem[j].dDistanceToCenter = w;
			VectorOfDistAndElem[j].iElementNumber = index;
			i++;
			j--;
		}
		else // when i >= j  we return j as the structure division point
			return j;
	}
}

//Overloaded partition function
int Partition(std::vector<AreaAndElemNumber> &VectorOfAreaAndElem, int p, int r) // we divide structure in 2 parts, in the first one all areas are less or equal to x, in the second one greater or equal to x
{
	double x = VectorOfAreaAndElem[p].dArea; // we take x
	int i = p, j = r, index;  // i, j - indexes in structure
	double w;
	while (true) // infinite loop - we leave it only by returning j
	{
		while (VectorOfAreaAndElem[j].dArea > x) // when elements are greater than x 
			j--;
		while (VectorOfAreaAndElem[i].dArea < x) // when elements are less than x
			i++;
		if (i < j) // we swap places when i < j
		{
			w = VectorOfAreaAndElem[i].dArea;
			index = VectorOfAreaAndElem[i].iElementNumber;

			VectorOfAreaAndElem[i].dArea = VectorOfAreaAndElem[j].dArea;
			VectorOfAreaAndElem[i].iElementNumber = VectorOfAreaAndElem[j].iElementNumber;

			VectorOfAreaAndElem[j].dArea = w;
			VectorOfAreaAndElem[j].iElementNumber = index;
			i++;
			j--;
		}
		else // when i >= j  we return j as the structure division point
			return j;
	}
}

//Overloaded quicksort algorithm
void Quicksort( std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r )
{
	int q;
	if (p < r)
	{
		q = Partition(VectorOfAngleAndElem, p, r); // we divide structure into 2 parts; q means separating point
		Quicksort(VectorOfAngleAndElem, p, q); // calling quicksort recursively for the first part of the structure
		Quicksort(VectorOfAngleAndElem, q + 1, r); // calling quicksort recursively for the second part of the structure
	}
}

//Overloaded quicksort algorithm
void Quicksort( std::vector<DistanceAndElemNumber> &VectorOfDistAndElem, int p, int r )
{
	int q;
	if (p < r)
	{
		q = Partition(VectorOfDistAndElem, p, r); // we divide structure into 2 parts; q means separating point
		Quicksort(VectorOfDistAndElem, p, q); // calling quicksort recursively for the first part of the structure
		Quicksort(VectorOfDistAndElem, q + 1, r); // calling quicksort recursively for the second part of the structure
	}
}

//Overloaded quicksort algorithm
void Quicksort(std::vector<AreaAndElemNumber> &VectorOfAreaAndElem, int p, int r)
{
	int q;
	if (p < r)
	{
		q = Partition(VectorOfAreaAndElem, p, r); // we divide structure into 2 parts; q means separating point
		Quicksort(VectorOfAreaAndElem, p, q); // calling quicksort recursively for the first part of the structure
		Quicksort(VectorOfAreaAndElem, q + 1, r); // calling quicksort recursively for the second part of the structure
	}
}

std::vector<Point2f> RemoveContoursFarestFromCenter( std::vector<DistanceAndElemNumber> VectorOfDistAndElem, std::vector<Point2f> mc2,
													 std::vector< vector <Point> > contours_reduced_2, std::vector< vector <Point> > &contours_reduced_3 )
{
	size_t iNumberOfPoints = VectorOfDistAndElem.size();
	std::vector<Point2f> mc3;
	int ElementNr;
	Point2f TempPoint;
	
	if ( iNumberOfPoints >= 4 )
	{
		for (size_t i = 0; i < 4; i++)
		{
			ElementNr = VectorOfDistAndElem[i].iElementNumber;
			TempPoint.x = mc2[ElementNr].x;
			TempPoint.y = mc2[ElementNr].y;
			mc3.push_back(TempPoint);

			vector<Point> contour = contours_reduced_2[ElementNr];
			contours_reduced_3.push_back(contour);
		}
	}

	return mc3;

}

std::vector<Point2f> EnumerateVerticiesClockwise( std::vector<AngleAndElemNumber> VectorOfAngleAndElem, std::vector<Point2f> mc3 )
{
	size_t iNumberOfPoints = VectorOfAngleAndElem.size();
	std::vector<Point2f> mc2_clockwise;
	Point2f TempPoint;
	int ElementNr;

	for (size_t i = 0; i < iNumberOfPoints; i++)
	{
		ElementNr = VectorOfAngleAndElem[i].iElementNumber;
		TempPoint.x = mc3[ElementNr].x;
		TempPoint.y = mc3[ElementNr].y;
		mc2_clockwise.push_back( TempPoint );
	}

	return mc2_clockwise;
}

std::vector< vector <Point> > RemoveSmallAndBigContours(std::vector< vector <Point> > contours, double dLowerAreaThreshold, double dHigherAreaThreshold)
{
	double area0;
	vector< vector <Point> > contours_reduced;

	for (size_t i = 0; i < contours.size(); i++)
	{
		vector<Point> contour = contours[i];
		area0 = contourArea(contour);
		if (area0 >= dLowerAreaThreshold && area0 <= dHigherAreaThreshold)
		{
			contours_reduced.push_back(contour);
		}
	}

	return contours_reduced;
}

std::vector< vector <Point> > RemoveObjEnclosingCircleLessThen( std::vector< vector <Point> > contours_reduced, double MinimumRadius )
{
	vector< vector <Point> > contours_reduced_2;
	Point2f center;
	float radius = 0.0;

	if ( contours_reduced.size() > 4 )
	{
		for (size_t i = 0; i < contours_reduced.size(); i++)
		{
			vector<Point> contour = contours_reduced[i];
			minEnclosingCircle(contour, center, radius);
			if ( radius > MinimumRadius )
			{
				contours_reduced_2.push_back(contour);
			}
		}
		return contours_reduced_2;
	}
	else
	{
		return contours_reduced;
	}	
}

std::vector<Point2f> GetContoursMassCenters(std::vector< vector <Point> > contours_reduced_2)
{
	///  Get the moments, from OpenCV documentation
	vector<Moments> mu(contours_reduced_2.size());
	vector<Point2f> mc2(contours_reduced_2.size());

	for (int i = 0; i < contours_reduced_2.size(); i++)
	{
		mu[i] = moments(contours_reduced_2[i], false);
	}

	///  Get the mass centers, from OpenCV documentation
	for (int i = 0; i < contours_reduced_2.size(); i++)
	{
		mc2[i] = Point2f(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00);
	}

	return mc2;
}


std::vector< vector <Point> > ApproximateContours( std::vector< vector <Point> > contours_reduced_2, double precision )
{
	vector< vector <Point> > contours_approxed;
	vector<Point> contour_approx;

	for (size_t i = 0; i < contours_reduced_2.size(); i++)
	{
		vector<Point> contour = contours_reduced_2[i];
		approxPolyDP(contour, contour_approx, precision, true);
		contours_approxed.push_back(contour_approx);
	}

	return contours_approxed;
}

double GetAverageNrOfPointsForContoursSet( std::vector< vector <Point> > contours_approxed )
{
	double dAvarageNumberOfPoints = 0.0;
	for (size_t i = 0; i < contours_approxed.size(); i++)
	{
		dAvarageNumberOfPoints = dAvarageNumberOfPoints + (double)contours_approxed[i].size();
	}

	dAvarageNumberOfPoints = dAvarageNumberOfPoints / (double)contours_approxed.size();
	return dAvarageNumberOfPoints;
}

std::vector<Point2f> GetTargetPointsDependingOnMarkerType( std::vector< vector <Point> > contours_approxed )
{
	Point2f Point1, Point2, Point3, Point4;
	std::vector<Point2f> TargetPoints;
	double dAvarageNumberOfPoints = GetAverageNrOfPointsForContoursSet( contours_approxed );
	if ( dAvarageNumberOfPoints < 5.0 )
	{
		Point1.x = 445; Point1.y = 0;
		Point2.x = 890; Point2.y = 315;
		Point3.x = 445; Point3.y = 630;
		Point4.x = 0; Point4.y = 315;
	}
	else
	{
		Point1.x = 890; Point1.y = 0;
		Point2.x = 890; Point2.y = 630;
		Point3.x = 0; Point3.y = 630;
		Point4.x = 0; Point4.y = 0;
	}

	TargetPoints.push_back(Point1);
	TargetPoints.push_back(Point2);
	TargetPoints.push_back(Point3);
	TargetPoints.push_back(Point4);

	return TargetPoints;
}

std::vector<AreaAndElemNumber> CalculateContoursArea( std::vector< vector <Point> > contours_jelly_reduced )
{
	std::vector<AreaAndElemNumber> VectorOfAreaAndElem;
	AreaAndElemNumber TempStruct;

	for (int i = 0; i < contours_jelly_reduced.size(); i++)
	{
		vector<Point> contour = contours_jelly_reduced[i];
		TempStruct.iElementNumber = i;
		TempStruct.dArea = contourArea( contour );

		VectorOfAreaAndElem.push_back(TempStruct);
	}

	return VectorOfAreaAndElem;
}

double Median( std::vector<AreaAndElemNumber> VectorOfAreaAndElem )
{
	double median = 0.0;
	int MiddleElementIndex = 0; //remember vector starts from 0 index, modified median algorithm needed
	size_t n = VectorOfAreaAndElem.size();
	if (n > 1)
	{
		if (n % 2 == 1)
		{
			MiddleElementIndex = (int)(n - 1) / 2;
			median = VectorOfAreaAndElem[MiddleElementIndex].dArea;
		}
		else
		{
			MiddleElementIndex = (int)n / 2;
			median = ( VectorOfAreaAndElem[MiddleElementIndex].dArea + VectorOfAreaAndElem[MiddleElementIndex - 1].dArea ) / 2;
		}
	}
	else if (n == 1)
	{
		median = VectorOfAreaAndElem[0].dArea;
	}

	return median;

}

//overloaded median function
int Median( std::vector<int> color_jellies )
{
	int median = 0, MiddleElementIndex = 0;  //remember vector starts from 0 index, modified median algorithm needed
	size_t n = color_jellies.size();

	if (n > 1)
	{
		if (n % 2 == 1)
		{
			MiddleElementIndex = (int)(n - 1) / 2;
			median = color_jellies[ MiddleElementIndex ];
		}
		else
		{
			MiddleElementIndex = (int)n / 2;
			median = (int)round( ( color_jellies[MiddleElementIndex] + color_jellies[MiddleElementIndex - 1] ) / 2 );
		}
	}
	else if (n == 1)
	{
		median = color_jellies[0];
	}

	return median;

}

int CountJellies( std::vector<AreaAndElemNumber> VectorOfAreaAndElem )
{
	int iNumberOfJellies = 0;
	double dMedian = Median( VectorOfAreaAndElem );
	if (dMedian > 5000)
	{
		dMedian = 3480; // experimentaly matched
	}
	if (VectorOfAreaAndElem.size() > 0)
	{
		for (size_t i = 0; i < VectorOfAreaAndElem.size(); i++)
		{
			iNumberOfJellies = iNumberOfJellies + (int)round( VectorOfAreaAndElem[i].dArea / dMedian );
		}
	}
	return iNumberOfJellies;
}

int CountColorJelly( cv::Mat warp_hsv, int color )
{
	Mat StructuringElement, dest_hsv, dilate_output, erode_output, light_red_secondary;
	int dilation_erosion_size = 3, iNumberOfJellies = 0;
	vector< vector <Point> > contours_jelly, contours_jelly_reduced;
	vector<AreaAndElemNumber> VectorOfAreaAndElem;
	Scalar LowerColorRange, UpperColorRange, LowerLightRedSecondary, HigherLightRedSecondary;

	// from OpenCV documentaion
	StructuringElement = getStructuringElement( MORPH_ELLIPSE,
		Size(2 * dilation_erosion_size + 1, 2 * dilation_erosion_size + 1),
		Point(dilation_erosion_size, dilation_erosion_size));

	switch (color)
	{
		case WHITE:
		{
			LowerColorRange = Scalar(11, 50, 103);
			UpperColorRange = Scalar(24, 150, 206);

			break;
		}
		case YELLOW:
		{
			LowerColorRange = Scalar(15, 163, 107);
			UpperColorRange = Scalar(47, 243, 205);
			break;
		}
		case ORANGE:
		{
			LowerColorRange = Scalar(6, 116, 82);
			UpperColorRange = Scalar(15, 255, 255);
			break;
		}
		case LIGHT_RED:
		{
			//Light red jelly in hsv consists of 2 colors (turquoise and green) - 2 ranges necessery
			LowerColorRange = Scalar(0, 165, 111);
			UpperColorRange = Scalar(5, 255, 255);

			LowerLightRedSecondary = Scalar(178, 152, 115);
			HigherLightRedSecondary = Scalar(180, 255, 255);
			break;
		}
		case DARK_RED:
		{
			LowerColorRange = Scalar(165, 115, 60);
			UpperColorRange = Scalar(179, 240, 132);
			break;
		}
		case GREEN:
		{
			LowerColorRange = Scalar(8, 77, 0);
			UpperColorRange = Scalar(56, 255, 124);
			break;
		}
	}
	
	inRange( warp_hsv, LowerColorRange, UpperColorRange, dest_hsv );
	//special case for light red jelly
	if ( color == LIGHT_RED )
	{
		inRange( warp_hsv, LowerLightRedSecondary, HigherLightRedSecondary, light_red_secondary );
		dest_hsv = dest_hsv + light_red_secondary;
	}
	//opening operation from OpenCV documentation
	erode( dest_hsv, erode_output, StructuringElement, cv::Point(-1, -1), 1 );
	dilate(erode_output, dilate_output, StructuringElement, cv::Point(-1, -1), 1);
	findContours(dilate_output, contours_jelly, RETR_TREE, CHAIN_APPROX_SIMPLE);
	contours_jelly_reduced = RemoveSmallAndBigContours( contours_jelly, 600, 25000 );

	if (contours_jelly_reduced.size() > 0)
	{
		VectorOfAreaAndElem = CalculateContoursArea(contours_jelly_reduced);
		Quicksort(VectorOfAreaAndElem, 0, (int)VectorOfAreaAndElem.size() - 1);
		iNumberOfJellies = CountJellies(VectorOfAreaAndElem);
	}

	contours_jelly.clear();
	contours_jelly_reduced.clear();
	VectorOfAreaAndElem.clear();

	return iNumberOfJellies;
}

int main(int, char)
{

	//reading text file - from tutorial: https://www.youtube.com/watch?v=h2Taf16gQDI
	vector<string> nazwy_obrazkow;
	vector<int> nr_sceny;
	int scena_int;
	string pojedyncza_linia, scena_str;
	fstream plik;
	plik.open("nazwy_zdjec/nazwy_zdjec.txt", ios::in);
	if ( plik.good() == false )
	{
		cout << "Podany plik nie istnieje lub podales bledna sciezke!" << endl;
		exit(0);
	}
	while ( getline( plik, pojedyncza_linia ) )
	{	
		nazwy_obrazkow.push_back( pojedyncza_linia );
		scena_str = pojedyncza_linia.substr(6, 3);
		scena_int = stoi(scena_str);
		nr_sceny.push_back(scena_int);
	}
	plik.close();


	//loading images to the vector<Mat> (with resizing) 
	vector <Mat> frames, frames_hsv;
	Mat temp_frame, temp_frame_grey, temp_frame_hsv;


	for ( size_t i = 0; i < nazwy_obrazkow.size(); i++ )
	{
		temp_frame = imread("zdjecia/" + nazwy_obrazkow[i], CV_LOAD_IMAGE_COLOR);
		cvtColor(temp_frame, temp_frame_grey, CV_BGR2GRAY);
		cvtColor(temp_frame, temp_frame_hsv, CV_BGR2HSV);
		Mat resized, resized2;

		if (temp_frame.rows <= temp_frame.cols)
		{
			resized.create(1170, 2080, temp_frame_grey.type());
			resized2.create(1170, 2080, temp_frame_hsv.type());
		}
		else
		{
			resized.create(2080, 1170, temp_frame_grey.type());
			resized2.create(2080, 1170, temp_frame_hsv.type());
		}
		
		cv::resize(temp_frame_grey, resized, resized.size());
		cv::resize(temp_frame_hsv, resized2, resized2.size());
		frames.push_back(resized);
		frames_hsv.push_back(resized2);
	}

	vector<int> white_jellies, yellow_jellies, orange_jellies, light_red_jellies, dark_red_jellies, green_jellies;
	fstream output_file;
	output_file.open("wyniki/Augustyn_Tomasz.txt", ios::out);

	for (size_t iCounter = 0; iCounter < frames.size(); iCounter++)

	{
		Mat threshold_img, StructuringElement, dilation_dst, warp_hsv, transform_matrix, erode_dst;
		int dilation_erosion_size = 3;
		vector< vector <Point> > contours, contours_reduced, contours_reduced_2, contours_reduced_3, contours_approxed;
		float radius = 0.0;
		double dLowerAreaThreshold = 0;
		vector<Point2f> mc2, mc3, mc3_related, mc3_clockwise;
		Point2f oCenterPoint;
		vector<AngleAndElemNumber> VectorOfAngleAndElem;
		vector<DistanceAndElemNumber> VectorOfDistAndElem;
		vector<Point2f> TargetPoints;


		for (int i = 0; i <= 4; i++)
		{
			threshold(frames[iCounter], threshold_img, 36 - i, 255, THRESH_BINARY);

			// from OpenCV documentaion
			StructuringElement = getStructuringElement(MORPH_ELLIPSE,
				Size(2 * dilation_erosion_size + 1, 2 * dilation_erosion_size + 1),
				Point(dilation_erosion_size, dilation_erosion_size));

			dilate(threshold_img, dilation_dst, StructuringElement, cv::Point(-1, -1), 3);
			erode(dilation_dst, erode_dst, StructuringElement, cv::Point(-1, -1), 2);

			findContours(erode_dst, contours, RETR_TREE, CHAIN_APPROX_SIMPLE);

			if (i == 0)
				dLowerAreaThreshold = 750;
			else
				dLowerAreaThreshold = 400;

			contours_reduced = RemoveSmallAndBigContours(contours, dLowerAreaThreshold, 18000);
			if (contours_reduced.size() <= 4)
			{
				break;
			}
		}

		contours_reduced_2 = RemoveObjEnclosingCircleLessThen(contours_reduced, 28);
		mc2 = GetContoursMassCenters(contours_reduced_2);

		if (contours_reduced_2.size() > 4)
		{
			VectorOfDistAndElem = CalculateDistanceToCenter(mc2, frames[iCounter]);
			Quicksort(VectorOfDistAndElem, 0, (int)VectorOfDistAndElem.size() - 1);
			mc3 = RemoveContoursFarestFromCenter(VectorOfDistAndElem, mc2, contours_reduced_2, contours_reduced_3);
		}
		else
		{
			contours_reduced_3 = contours_reduced_2;
			mc3 = mc2;
		}

		oCenterPoint = CalculateCenterPointOfPoints(mc3);
		mc3_related = RelateVectorOfPointsToTheCenterPoint(mc3, oCenterPoint);
		VectorOfAngleAndElem = ConvertToPolarCoordinates(mc3_related);
		Quicksort(VectorOfAngleAndElem, 0, (int)VectorOfAngleAndElem.size() - 1);
		mc3_clockwise = EnumerateVerticiesClockwise(VectorOfAngleAndElem, mc3);

		int iWhiteJellies = 0, iYellowJellies = 0, iOrangeJellies = 0, iLightRedJellies = 0, iDarkRedJellies = 0, iGreenJellies = 0;
		int MiddleWhite = 0, MiddleYellow = 0, MiddleOrange = 0, MiddleLightRed = 0, MiddleDarkRed = 0, MiddleGreen = 0;


		if (mc3_clockwise.size() == 4 && contours_reduced_3.size() == 4)
		{
			contours_approxed = ApproximateContours(contours_reduced_3, 8);
			TargetPoints = GetTargetPointsDependingOnMarkerType(contours_approxed);
			transform_matrix = getPerspectiveTransform(mc3_clockwise, TargetPoints);
			warpPerspective(frames_hsv[iCounter], warp_hsv, transform_matrix, Size(890, 630));

			iWhiteJellies = CountColorJelly(warp_hsv, WHITE);
			iYellowJellies = CountColorJelly(warp_hsv, YELLOW);
			iOrangeJellies = CountColorJelly(warp_hsv, ORANGE);
			iLightRedJellies = CountColorJelly(warp_hsv, LIGHT_RED);
			iDarkRedJellies = CountColorJelly(warp_hsv, DARK_RED);
			iGreenJellies = CountColorJelly(warp_hsv, GREEN);

			if (iCounter != 0)
			{
				if (nr_sceny[iCounter] == nr_sceny[iCounter - 1])
				{
					white_jellies.push_back(iWhiteJellies);
					yellow_jellies.push_back(iYellowJellies);
					orange_jellies.push_back(iOrangeJellies);
					light_red_jellies.push_back(iLightRedJellies);
					dark_red_jellies.push_back(iDarkRedJellies);
					green_jellies.push_back(iGreenJellies);
				}
				else
				{
					std::sort(white_jellies.begin(), white_jellies.end());
					std::sort(yellow_jellies.begin(), yellow_jellies.end());
					std::sort(orange_jellies.begin(), orange_jellies.end());
					std::sort(light_red_jellies.begin(), light_red_jellies.end());
					std::sort(dark_red_jellies.begin(), dark_red_jellies.end());
					std::sort(green_jellies.begin(), green_jellies.end());

					MiddleWhite = Median(white_jellies);
					MiddleYellow = Median(yellow_jellies);
					MiddleOrange = Median(orange_jellies);
					MiddleLightRed = Median(light_red_jellies);
					MiddleDarkRed = Median(dark_red_jellies);
					MiddleGreen = Median(green_jellies);

					output_file << MiddleDarkRed << ", " << MiddleLightRed << ", " << MiddleGreen << ", " << MiddleOrange << ", " << MiddleWhite << ", " << MiddleYellow << endl;

					white_jellies.clear();
					yellow_jellies.clear();
					orange_jellies.clear();
					light_red_jellies.clear();
					dark_red_jellies.clear();
					green_jellies.clear();

					white_jellies.push_back(iWhiteJellies);
					yellow_jellies.push_back(iYellowJellies);
					orange_jellies.push_back(iOrangeJellies);
					light_red_jellies.push_back(iLightRedJellies);
					dark_red_jellies.push_back(iDarkRedJellies);
					green_jellies.push_back(iGreenJellies);

				}

				if ( iCounter == (frames.size() - 1) )
				{
					std::sort(white_jellies.begin(), white_jellies.end());
					std::sort(yellow_jellies.begin(), yellow_jellies.end());
					std::sort(orange_jellies.begin(), orange_jellies.end());
					std::sort(light_red_jellies.begin(), light_red_jellies.end());
					std::sort(dark_red_jellies.begin(), dark_red_jellies.end());
					std::sort(green_jellies.begin(), green_jellies.end());

					MiddleWhite = Median(white_jellies);
					MiddleYellow = Median(yellow_jellies);
					MiddleOrange = Median(orange_jellies);
					MiddleLightRed = Median(light_red_jellies);
					MiddleDarkRed = Median(dark_red_jellies);
					MiddleGreen = Median(green_jellies);

					output_file << MiddleDarkRed << ", " << MiddleLightRed << ", " << MiddleGreen << ", " << MiddleOrange << ", " << MiddleWhite << ", " << MiddleYellow << endl;
				}

			}
			else
			{

				white_jellies.push_back( iWhiteJellies );
				yellow_jellies.push_back( iYellowJellies );
				orange_jellies.push_back( iOrangeJellies );
				light_red_jellies.push_back( iLightRedJellies );
				dark_red_jellies.push_back( iDarkRedJellies );
				green_jellies.push_back( iGreenJellies );

				if ( frames.size() == 1 )
				{
					output_file << MiddleDarkRed << ", " << MiddleLightRed << ", " << MiddleGreen << ", " << MiddleOrange << ", " << MiddleWhite << ", " << MiddleYellow << endl;
				}
			}

		}
		else if (iCounter == (frames.size() - 1))
		{
			std::sort(white_jellies.begin(), white_jellies.end());
			std::sort(yellow_jellies.begin(), yellow_jellies.end());
			std::sort(orange_jellies.begin(), orange_jellies.end());
			std::sort(light_red_jellies.begin(), light_red_jellies.end());
			std::sort(dark_red_jellies.begin(), dark_red_jellies.end());
			std::sort(green_jellies.begin(), green_jellies.end());

			MiddleWhite = Median(white_jellies);
			MiddleYellow = Median(yellow_jellies);
			MiddleOrange = Median(orange_jellies);
			MiddleLightRed = Median(light_red_jellies);
			MiddleDarkRed = Median(dark_red_jellies);
			MiddleGreen = Median(green_jellies);

			output_file << MiddleDarkRed << ", " << MiddleLightRed << ", " << MiddleGreen << ", " << MiddleOrange << ", " << MiddleWhite << ", " << MiddleYellow << endl;
		}

	}

	output_file.close();
	cout << "Zapis wynikow do pliku zakonczony."<<endl;

	white_jellies.clear();
	yellow_jellies.clear();
	orange_jellies.clear();
	light_red_jellies.clear();
	dark_red_jellies.clear();
	green_jellies.clear();

	frames.clear();
	frames_hsv.clear();
	nazwy_obrazkow.clear();
	nr_sceny.clear();

	waitKey();
	return 0;
}