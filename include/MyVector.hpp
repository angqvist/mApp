

class MyVector
{
public:
  myVector();
  myVector(int);
  friend MyVector operator+(MyVector,MyVector);
  friend MyVector operator+=(MyVector,MyVector);
  friend MyVector operator-(MyVector,MyVector);
  friend MyVector operator-=(MyVector,MyVector);
  friend MyVector operator*(MyVector,MyVector);	 
  double& operator[](int);
  const double& operator[](int) const;

private:
  std::vector<double> theVector;
  int length;
  
