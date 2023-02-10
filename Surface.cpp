#ifndef SURFACE
#define SURFACE

#include "Position.cpp"
#include "math.h"
#include <limits>

class Surface {
    protected:
        Position p0;
        std::string type;

    public:
        Surface(Position p0):p0(p0) {
            type = "surface";
        }

        Position getLocation() {
            return p0;
        }

        virtual bool senseOfPosition(Position p) = 0;

        virtual double distanceToSurface(Position location, Direction v) = 0;
};

class Plane: public Surface {
    protected:
        Direction norm;

    public:
        Plane(Position p0, Direction norm):Surface(p0),norm(norm) {
            type = "plane";
        }
        
        Direction getNormalVec() {
            return norm;
        }

        // returns true if position is on the normal vector side of the plane (sense of 1) 
        // and false if the position is on the opposite side of the plane (sense of 0)
        bool senseOfPosition(Position p) override{
            double solvedEquation = norm.getI() * (p.x - p0.x) + norm.getJ() * (p.y - p0.y) + norm.getK() * (p.z - p0.z);
            return solvedEquation > 0;
        }

        double distanceToSurface(Position location, Direction v) override{
            double d = (norm.getI() * (location.getX() - p0.getX()) + norm.getJ() * (location.getY() - p0.getY()) + norm.getK() * (location.getZ() - p0.getZ())) / (norm.getI() * v.getI() + norm.getJ() * v.getJ() + norm.getK() * v.getK());
            if (d < 0) {
                return std::numeric_limits<double>::max();
            }
            return d;
        }
};

class Sphere: public Surface {
    protected:
        double radius;

    public:
        Sphere(Position p0, double radius):Surface(p0),radius(radius) {
            type = "sphere";
        }

        // return true if position is outside sphere (sense of 1) and false if position is inside sphere (sense of 0)
        bool senseOfPosition(Position p) override{
            double distFromCenter = sqrt(pow(p.x-p0.x,2) + pow(p.y-p0.y,2) + pow(p.z-p0.z,2));
            return distFromCenter > radius;
        }

        double distanceToSurface(Position location, Direction v) override{
            double a = pow(v.getI(),2) + pow(v.getJ(),2) + pow(v.getK(),2);
            double b = 2 * (v.getI() * (location.getX() - p0.getX()) + v.getJ() * (location.getY() - p0.getY()) + v.getK() * (location.getZ() - p0.getZ()));
            double c = pow(location.getX() - p0.getX(), 2) + pow(location.getY() - p0.getY(), 2) + pow(location.getZ() - p0.getZ(), 2) - pow(radius,2);

            double discriminant = pow(b,2) - 4 * a * c;

            if(discriminant < 0) {
                return std::numeric_limits<double>::max();
            }
            else if (-1 * b - sqrt(discriminant) > 0){
                return (-1 * b - sqrt(discriminant)) / (2 * a);
            }
            else if (-1 * b + sqrt(discriminant) > 0){
                return (-1 * b + sqrt(discriminant)) / (2 * a);
            }
            return std::numeric_limits<double>::max();
        }
};

#endif