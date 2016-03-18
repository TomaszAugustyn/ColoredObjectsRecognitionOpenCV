#include "opencv2/opencv.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>


using namespace cv;
using namespace std;

#define M_PI           3.14159265358979323846  /* pi */

//declarations
struct AngleAndElemNumber;
struct DistanceAndElemNumber;
Point2f CalculateCenterPointOfPoints(std::vector<Point2f> mc2);
std::vector<Point2f> RelateVectorOfPointsToTheCenterPoint(std::vector<Point2f> mc2, Point2f oCenterPoint);
std::vector<AngleAndElemNumber> ConvertToPolarCoordinates(std::vector<Point2f> mc2_related);
std::vector<DistanceAndElemNumber> CalculateDistanceToCenter(std::vector<Point2f> mc2, cv::Mat frame);
//int Partition(std::vector< vector <Point> > &contours_reduced_2, int p, int r);
int Partition(std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r);
int Partition(std::vector<DistanceAndElemNumber> &VectorOfDistAndElem, int p, int r);
//void Quicksort(std::vector< vector <Point> > &contours_reduced_2, int p, int r);
void Quicksort(std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r);
void Quicksort(std::vector<DistanceAndElemNumber> &VectorOfDistAndElem, int p, int r);
std::vector<Point2f> RemoveContoursFarestFromCenter(std::vector<DistanceAndElemNumber> VectorOfDistAndElem, std::vector<Point2f> mc2, std::vector< vector <Point> > contours_reduced_2, std::vector< vector <Point> > &contours_reduced_3);
std::vector<Point2f> EnumerateVerticiesClockwise(std::vector<AngleAndElemNumber> VectorOfAngleAndElem, std::vector<Point2f> mc2);
std::vector< vector <Point> > RemoveSmallAndBigContours(std::vector< vector <Point> > contours, double dLowerAreaThreshold, double dHigherAreaThreshold);
std::vector< vector <Point> > RemoveObjEnclosingCircleLessThen(std::vector< vector <Point> > contours_reduced, double MinimumRadius);
//void RemoveTooComplicatedContours(std::vector< vector <Point> > &contours_reduced_2);
std::vector< vector <Point> > ApproximateContours(std::vector< vector <Point> > contours_reduced_2, double precision);
double GetAverageNrOfPointsForContoursSet(std::vector< vector <Point> > contours_approxed);
std::vector<Point2f> GetTargetPointsDependingOnMarkerType(std::vector< vector <Point> > contours_approxed);

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

Point2f CalculateCenterPointOfPoints(std::vector<Point2f> mc3)
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

std::vector<Point2f> RelateVectorOfPointsToTheCenterPoint(std::vector<Point2f> mc3, Point2f oCenterPoint)
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

std::vector<AngleAndElemNumber> ConvertToPolarCoordinates(std::vector<Point2f> mc3_related)
{
	int iNumberOfPoints = mc3_related.size();
	std::vector<AngleAndElemNumber> VectorOfAngleAndElem;
	AngleAndElemNumber TempStruct;
	float tempX = 0.0, tempY = 0.0;

	for (int i = 0; i < iNumberOfPoints; i++)
	{
		TempStruct.iElementNumber = i;
		tempX = mc3_related[i].x;
		tempY = mc3_related[i].y;

		if ((tempX > 0) && (tempY >= 0))
		{
			TempStruct.dPolarAngle = atan(tempY / tempX);
		}
		else if ((tempX > 0) && (tempY < 0))
		{
			TempStruct.dPolarAngle = atan(tempY / tempX) + 2 * M_PI;
		}
		else if (tempX < 0)
		{
			TempStruct.dPolarAngle = atan(tempY / tempX) + M_PI;
		}
		else if ((tempX == 0) && (tempY > 0))
		{
			TempStruct.dPolarAngle = M_PI / 2;
		}
		else if ((tempX == 0) && (tempY < 0))
		{
			TempStruct.dPolarAngle = 3 * M_PI / 2;
		}
		else
		{
			TempStruct.dPolarAngle = 0;
		}

		VectorOfAngleAndElem.push_back(TempStruct);
	}

	return VectorOfAngleAndElem;
}

std::vector<DistanceAndElemNumber> CalculateDistanceToCenter(std::vector<Point2f> mc2, cv::Mat frame)
{
	int iNumberOfPoints = mc2.size();
	std::vector<DistanceAndElemNumber> VectorOfDistAndElem;
	DistanceAndElemNumber TempStruct;
	double tempX = 0.0, tempY = 0.0;
	double CenterX = frame.cols / 2;
	double CenterY = frame.rows / 2;
	//float CenterX = 1040;
	//float CenterY = 585;
	double distance = 0.0;


	for (int i = 0; i < iNumberOfPoints; i++)
	{
		TempStruct.iElementNumber = i;
		tempX = mc2[i].x;
		tempY = mc2[i].y;
		TempStruct.dDistanceToCenter = sqrt(pow(CenterX - tempX, 2) + pow(CenterY - tempY, 2));

		VectorOfDistAndElem.push_back(TempStruct);
	}

	return VectorOfDistAndElem;
}

/*int Partition(std::vector< vector <Point> > &contours_reduced_2, int p, int r) // dzielimy wektor na dwie czesci, w pierwszej rozmiar wektora punktów jest mniejszy badz rowny x, w drugiej wiekszy lub rowny od x
{
int x = contours_reduced_2[p].size(); // obieramy x
int i = p, j = r; // i, j - indeksy wektora konturow
vector <Point> w; // pojedynczy kontur (wektor punktow)
while (true) // petla nieskonczona - wychodzimy z niej tylko przez return j
{
while (contours_reduced_2[j].size() > x) // dopoki elementy sa wieksze od x
j--;
while (contours_reduced_2[i].size() < x) // dopoki elementy sa mniejsze od x
i++;
if (i < j) // zamieniamy miejscami gdy i < j
{
w = contours_reduced_2[i];
contours_reduced_2[i] = contours_reduced_2[j];
contours_reduced_2[j] = w;
i++;
j--;
}
else // gdy i >= j zwracamy j jako punkt podzialu tablicy
return j;
}
}*/

//Overloaded partition function
int Partition(std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r) // dzielimy strukturê na dwie czesci, w pierwszej wszystkie rozmiary katow sa mniejsze badz rowne x, w drugiej wieksze lub rowne od x
{
	double x = VectorOfAngleAndElem[p].dPolarAngle; // obieramy x
	int i = p, j = r, index; // i, j - indeksy w strukturze
	double w;
	while (true) // petla nieskonczona - wychodzimy z niej tylko przez return j
	{
		while (VectorOfAngleAndElem[j].dPolarAngle > x) // dopoki elementy sa wieksze od x
			j--;
		while (VectorOfAngleAndElem[i].dPolarAngle < x) // dopoki elementy sa mniejsze od x
			i++;
		if (i < j) // zamieniamy miejscami gdy i < j
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
		else // gdy i >= j zwracamy j jako punkt podzialu tablicy
			return j;
	}
}

//Overloaded partition function
int Partition(std::vector<DistanceAndElemNumber> &VectorOfDistAndElem, int p, int r) // dzielimy strukturê na dwie czesci, w pierwszej wszystkie rozmiary katow sa mniejsze badz rowne x, w drugiej wieksze lub rowne od x
{
	double x = VectorOfDistAndElem[p].dDistanceToCenter; // obieramy x
	int i = p, j = r, index; // i, j - indeksy w strukturze
	double w;
	while (true) // petla nieskonczona - wychodzimy z niej tylko przez return j
	{
		while (VectorOfDistAndElem[j].dDistanceToCenter > x) // dopoki elementy sa wieksze od x
			j--;
		while (VectorOfDistAndElem[i].dDistanceToCenter < x) // dopoki elementy sa mniejsze od x
			i++;
		if (i < j) // zamieniamy miejscami gdy i < j
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
		else // gdy i >= j zwracamy j jako punkt podzialu tablicy
			return j;
	}
}


/*void Quicksort(std::vector< vector <Point> > &contours_reduced_2, int p, int r) // sortowanie szybkie
{
int q;
if (p < r)
{
q = Partition(contours_reduced_2, p, r); // dzielimy wektor konturow na dwie czesci; q oznacza punkt podzialu
Quicksort(contours_reduced_2, p, q); // wywolujemy rekurencyjnie quicksort dla pierwszej czesci wektora
Quicksort(contours_reduced_2, q + 1, r); // wywolujemy rekurencyjnie quicksort dla drugiej czesci wektora
}
}*/

//Overloaded quicksort algorithm
void Quicksort(std::vector<AngleAndElemNumber> &VectorOfAngleAndElem, int p, int r) // sortowanie szybkie
{
	int q;
	if (p < r)
	{
		q = Partition(VectorOfAngleAndElem, p, r); // dzielimy strukture na dwie czesci; q oznacza punkt podzialu
		Quicksort(VectorOfAngleAndElem, p, q); // wywolujemy rekurencyjnie quicksort dla pierwszej czesci struktury
		Quicksort(VectorOfAngleAndElem, q + 1, r); // wywolujemy rekurencyjnie quicksort dla drugiej czesci struktury
	}
}

//Overloaded quicksort algorithm
void Quicksort(std::vector<DistanceAndElemNumber> &VectorOfDistAndElem, int p, int r) // sortowanie szybkie
{
	int q;
	if (p < r)
	{
		q = Partition(VectorOfDistAndElem, p, r); // dzielimy strukture na dwie czesci; q oznacza punkt podzialu
		Quicksort(VectorOfDistAndElem, p, q); // wywolujemy rekurencyjnie quicksort dla pierwszej czesci struktury
		Quicksort(VectorOfDistAndElem, q + 1, r); // wywolujemy rekurencyjnie quicksort dla drugiej czesci struktury
	}
}

std::vector<Point2f> RemoveContoursFarestFromCenter(std::vector<DistanceAndElemNumber> VectorOfDistAndElem, std::vector<Point2f> mc2,
	std::vector< vector <Point> > contours_reduced_2, std::vector< vector <Point> > &contours_reduced_3)
{
	size_t iNumberOfPoints = VectorOfDistAndElem.size();
	std::vector<Point2f> mc3;
	int ElementNr;
	Point2f TempPoint;

	if (iNumberOfPoints >= 4)
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

std::vector<Point2f> EnumerateVerticiesClockwise(std::vector<AngleAndElemNumber> VectorOfAngleAndElem, std::vector<Point2f> mc3)
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
		mc2_clockwise.push_back(TempPoint);
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

std::vector< vector <Point> > RemoveObjEnclosingCircleLessThen(std::vector< vector <Point> > contours_reduced, double MinimumRadius)
{
	vector< vector <Point> > contours_reduced_2;
	Point2f center;
	float radius = 0.0;

	if (contours_reduced.size() > 4)
	{
		for (size_t i = 0; i < contours_reduced.size(); i++)
		{
			vector<Point> contour = contours_reduced[i];
			minEnclosingCircle(contour, center, radius);
			if (radius > MinimumRadius)
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

/*void RemoveTooComplicatedContours( std::vector< vector <Point> > &contours_reduced_2 )
{
if (contours_reduced_2.size() > 4)
{
Quicksort( contours_reduced_2, 0, (int)contours_reduced_2.size() - 1 );
for ( size_t i = contours_reduced_2.size() - 1; i >= 4; i-- )
{
contours_reduced_2.erase( contours_reduced_2.begin() + i );
}
}
}*/

std::vector< vector <Point> > ApproximateContours(std::vector< vector <Point> > contours_reduced_2, double precision)
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

double GetAverageNrOfPointsForContoursSet(std::vector< vector <Point> > contours_approxed)
{
	double dAvarageNumberOfPoints = 0.0;
	for (size_t i = 0; i < contours_approxed.size(); i++)
	{
		dAvarageNumberOfPoints = dAvarageNumberOfPoints + (double)contours_approxed[i].size();
	}

	dAvarageNumberOfPoints = dAvarageNumberOfPoints / (double)contours_approxed.size();
	return dAvarageNumberOfPoints;
}

std::vector<Point2f> GetTargetPointsDependingOnMarkerType(std::vector< vector <Point> > contours_approxed)
{
	Point2f p1, p2, p3, p4;
	std::vector<Point2f> TargetPoints;
	double dAvarageNumberOfPoints = GetAverageNrOfPointsForContoursSet(contours_approxed);
	if (dAvarageNumberOfPoints < 5.0)
	{
		p1.x = 445;
		p1.y = 0;

		p2.x = 890;
		p2.y = 315;

		p3.x = 445;
		p3.y = 630;

		p4.x = 0;
		p4.y = 315;
	}
	else
	{
		p1.x = 890;
		p1.y = 0;

		p2.x = 890;
		p2.y = 630;

		p3.x = 0;
		p3.y = 630;

		p4.x = 0;
		p4.y = 0;
	}

	TargetPoints.push_back(p1);
	TargetPoints.push_back(p2);
	TargetPoints.push_back(p3);
	TargetPoints.push_back(p4);

	return TargetPoints;
}

int main(int, char)
{
	//odczyt pliku tekstowego
	vector<string> nazwy_obrazkow;
	string pojedyncza_linia;
	fstream plik;
	plik.open("nazwy_zdjec/nazwy_zdjec.txt", ios::in);
	if (plik.good() == false)
	{
		cout << "Podany plik nie istnieje lub podales bledna sciezke!" << endl;
		exit(0);
	}
	while (getline(plik, pojedyncza_linia))
	{
		nazwy_obrazkow.push_back(pojedyncza_linia);
	}
	plik.close();


	//ladowanie obrazkow do wektora Mat i skalowanie ich
	vector <Mat> frames;
	Mat temp_frame;
	for (size_t i = 0; i < nazwy_obrazkow.size(); i++)
	{
		temp_frame = imread("zdjecia/" + nazwy_obrazkow[i], CV_LOAD_IMAGE_GRAYSCALE);
		Mat resized;

		if (temp_frame.rows <= temp_frame.cols)
		{
			resized.create(1170, 2080, temp_frame.type());

		}
		else
		{
			resized.create(2080, 1170, temp_frame.type());
		}

		cv::resize(temp_frame, resized, resized.size());
		frames.push_back(resized);
	}


	Mat threshold_img, element, dilation_dst, dest_img, transform_matrix, erode_dst, test;
	int dilation_size = 3;
	vector< vector <Point> > contours, contours_reduced, contours_reduced_2, contours_reduced_3, contours_approxed;
	vector <double> areas;
	double area0 = 0;
	Point2f center;
	float radius = 0.0;
	double dLowerAreaThreshold = 0;
	namedWindow("okno", 1);
	namedWindow("okno2", 1);

	for (int i = 0; i <= 4; i++)
	{
		threshold(frames[13], threshold_img, 36 - i, 255, THRESH_BINARY);

		element = getStructuringElement(MORPH_ELLIPSE,
			Size(2 * dilation_size + 1, 2 * dilation_size + 1),
			Point(dilation_size, dilation_size));

		dilate(threshold_img, dilation_dst, element, cv::Point(-1, -1), 3);
		erode(dilation_dst, erode_dst, element, cv::Point(-1, -1), 2);

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
	//	RemoveTooComplicatedContours( contours_reduced_2 );

	for (size_t i = 0; i < contours_reduced_2.size(); i++)
	{
		vector<Point> contour = contours_reduced_2[i];

		area0 = contourArea(contour);
		areas.push_back(area0);
	}

	/// Get the moments
	vector<Moments> mu(contours_reduced_2.size());
	for (int i = 0; i < contours_reduced_2.size(); i++)
	{
		mu[i] = moments(contours_reduced_2[i], false);
	}

	///  Get the mass centers:
	vector<Point2f> mc2(contours_reduced_2.size());
	vector<Point2f> mc3, mc3_related, mc3_clockwise;
	Point2f oCenterPoint;
	vector<AngleAndElemNumber> VectorOfAngleAndElem;
	vector<DistanceAndElemNumber> VectorOfDistAndElem;
	vector<Point2f> TargetPoints;

	for (int i = 0; i < contours_reduced_2.size(); i++)
	{
		mc2[i] = Point2f(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00);
	}

	if (contours_reduced_2.size() > 4)
	{
		VectorOfDistAndElem = CalculateDistanceToCenter(mc2, frames[13]);
		Quicksort(VectorOfDistAndElem, 0, (int)VectorOfDistAndElem.size() - 1);
		mc3 = RemoveContoursFarestFromCenter(VectorOfDistAndElem, mc2, contours_reduced_2, contours_reduced_3);
	}
	else
	{
		contours_reduced_3 = contours_reduced_2;
		mc3 = mc2;
	}

	Scalar color(255, 255, 255);
	for (int i = 1; i < contours_reduced_3.size(); i++)
	{
		//drawContours(dilation_dst, contours, i, color, CV_FILLED);
		drawContours(erode_dst, contours_reduced_3, i, color, 2, 8);
		Mat resized2;
		resized2.create(frames[13].rows / 2, frames[13].cols / 2, frames[7].type());
		cv::resize(erode_dst, test, resized2.size());

	}

	oCenterPoint = CalculateCenterPointOfPoints(mc3);
	mc3_related = RelateVectorOfPointsToTheCenterPoint(mc3, oCenterPoint);
	VectorOfAngleAndElem = ConvertToPolarCoordinates(mc3_related);
	Quicksort(VectorOfAngleAndElem, 0, (int)VectorOfAngleAndElem.size() - 1);
	mc3_clockwise = EnumerateVerticiesClockwise(VectorOfAngleAndElem, mc3);

	contours_approxed = ApproximateContours(contours_reduced_3, 7);
	TargetPoints = GetTargetPointsDependingOnMarkerType(contours_approxed);




	transform_matrix = getPerspectiveTransform(mc3_clockwise, TargetPoints);
	/*if (frames[13].rows <= frames[13].cols)
	{*/
	warpPerspective(frames[13], dest_img, transform_matrix, Size(890, 630));
	/*	}
	else
	{
	warpPerspective(frames[13], dest_img, transform_matrix, Size(500, 800));
	}*/
	imshow("okno2", test);
	imshow("okno", dest_img);

	waitKey();
	return 0;
}