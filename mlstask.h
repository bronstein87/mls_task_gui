#ifndef MLSTASK_H
#define MLSTASK_H

#include <QObject>
#include <QFile>
#include <QTextStream>
#include <QVector>
#include <QPointF>
#include <QBitArray>
#include <mathfunc.h>
#include <gms.h>
using namespace std;
using namespace BOKZMath;
// изначально передаются первые приближения
struct Results
{
    double foc = 0;
    double lambdaOY = 0;
    double lambdaOX = 0;
    double alphaTwoOY = 0;
    double alphaTwoOX = 0;
    double gammaOZ = 0;
    double gammaOY = 0;
    double gammaOX = 0;
};

struct ResultErrors
{
    double dfoc = 0;
    double dlambdaOY = 0;
    double dlambdaOX = 0;
    double dalphaTwoOY = 0;
    double dalphaTwoOX = 0;
    double dgammaOZ = 0;
    double dgammaOY = 0;
    double dgammaOX = 0;
    double mx = 0;
    double my = 0;
    double mxy = 0;
};

struct RotateAngles
{
    QVector <double> alphaRotates;
    QVector <double> phiRotates;
};


enum DERIVATIVES
{
    LAMBDA_OY,
    LAMBDA_OX,
    ALPHA_TWO_OY,
    ALPHA_TWO_OX,
    GAMMA_OZ,
    GAMMA_OY,
    GAMMA_OX,
    PHI,
    FOCUS,
    COUNT

};

struct DistorsioData
{
    QVector <double> x;
    QVector <double> y;
    QVector <double> dx;
    QVector <double> dy;
    QVector <double> dx_diff;
    QVector <double> dy_diff;
};


enum X_AXIS_DIRECTION
{
    UP,
    DOWN,
    LEFT,
    RIGHT,
    REVERSE_X,
    REVERSE_Y
};

class MLSTask : public QObject
{
    Q_OBJECT

    struct StandAngles
    {
        double lambdaOY = 0;
        double lambdaOX = 0;
        double alphaTwoOY = 0;
        double alphaTwoOX = 0;
        double gammaOZ =0;
        double gammaOY = 0;
        double gammaOX = 0;
    };

public:

    explicit MLSTask(QObject* parent = nullptr);
    void readModelData(const QString& filename, bool skipFirstRow);
    void readRealData(const QString& filename, bool skipFirstRow, bool reverse = false);
    void calculate(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors);
    void fitFocusByLines(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors);
    void setPixelSize(double _pixelSize) {pixelSize = _pixelSize;}
    void setFrameSize(quint32 x, quint32 y) {frameX = x; frameY = y;}
    void setMeasureTheshold(quint32 th){thershold = th;}
    void findDistorsio(int nPow);
    void saveDistorsio();
    void saveShifts(const QString& prefix);
    void includeDistorsio();
    QVector<QString> printTestTable(const QString& filename, bool dist, double focus);
    QList <double> getDistX() const {return xDistV;}
    QList <double> getDistY() const {return yDistV;}
    void setAxisDirection(X_AXIS_DIRECTION dir) {d = dir;}
    void clearAll(){frame.clear();xDistV.clear();yDistV.clear();distData.x.clear();distData.y.clear();distData.dx.clear();distData.dy.clear();}

private:
    void includeAxisDirection (double MStand[3][3], double modifyMStand[3][3], X_AXIS_DIRECTION d = X_AXIS_DIRECTION::UP);
    void calculatePrivate(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors, const RotateAngles& rotAngles, QVector<QPointF>& frame);
    void FindDistCft(int Npow, QVector <double>& x, QVector <double>& y, QVector <double>& dx, QVector <double>& dy);
    void calculateXY(double Mstand[3][3], double collimator[3], double focus, double& X, double& Y);
    void calculateMatrix(StandAngles ang, double alpha, double phi, double Mstand[3][3]);
    void calculateMatrixDynamicaly(StandAngles ang, double alpha, double phi, double Mstand[3][3]);
    double calculate11El(StandAngles ang, double alpha, double phi);
    double calculate21El(StandAngles ang, double alpha, double phi);
    double calculate31El(StandAngles ang, double alpha, double phi);
    double calculate12El(StandAngles ang, double alpha, double phi);
    double calculate22El(StandAngles ang, double alpha, double phi);
    double calculate32El(StandAngles ang, double alpha, double phi);
    double calculate13El(StandAngles ang, double alpha, double phi);
    double calculate23El(StandAngles ang, double alpha, double phi);
    double calculate33El(StandAngles ang, double alpha, double phi);
    int gaus_obr(int cnt_str, double mass[55][55], double M_obr[55][55]);

    double pixelSize;
    quint32 thershold = 2;
    quint32 frameX;
    quint32 frameY;
    QVector <QPointF> frame;
    QVector <QPointF> initFrame;
    RotateAngles rotAngles;
    DistorsioData distData;
    QList <double> xDistV;
    QList <double> yDistV;
    X_AXIS_DIRECTION d = UP;
    static constexpr const double deltaAngle = 1. / 60. / 60. * degreesToRad * 10;
    static constexpr const double deltaFocus = 0.0001;
    static constexpr const qint32 maxDer = 10000;
    static constexpr const qint32 maxParams = 55;
    static constexpr const qint32 maxIterations = 500;
    static constexpr const double maxError = 0.000005;
};

#endif // MLSTASK_H
