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
#include <functional>
using namespace std;
using namespace BOKZMath;

struct Catalog
{
    QVector <double> L;
    QVector <double> M;
    QVector <double> N;
    QVector <double> alpha;
    QVector <double> delta;
    quint32 countStars = 0;
};
struct OrientAngles
{
    double   al;
    double   dl;
    double   Az;
};

struct OldErrors
{
    double dAl;
    double dDlt;
    double dAz;
    double dFoc;
};

struct LineData
{
    QVector <double> x;
    QVector <double> y;
    QVector <qint32> vecSt;
    qint32 count = 0;
    void clear() {x.clear();y.clear(); vecSt.clear(); count = 0;}
};

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
    QVector <double> dxDiff;
    QVector <double> dyDiff;
};

struct ShiftData
{
    QVector <double> x;
    QVector <double> y;
    QVector <double> dx;
    QVector <double> dy;
    QVector <double> azimut;
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
static constexpr const qint32 countPointsTriangle = 3;
struct TrianglePattern
{

    TrianglePattern () : points(countPointsTriangle) {}

    //    TrianglePattern(TrianglePattern&& fr) : firstPoint(qMove(fr.firstPoint)),
    //        secondPoint(qMove(fr.secondPoint)), thirdPoint(qMove(fr.thirdPoint)),
    //        points(qMove(fr.points))
    //    {}

    //    TrianglePattern(TrianglePattern& fr) : firstPoint(fr.firstPoint),
    //        secondPoint(fr.secondPoint), thirdPoint(fr.thirdPoint),
    //        points(fr.points){}

    //    QPointF  firstPoint;
    //    QPointF  secondPoint;
    //    QPointF  thirdPoint;
    /* const*/ //QVector <std::reference_wrapper<QPointF>> points;
    QVector <QPointF> points;
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

    void readFullModelData(const QString& filename, bool skipFirstRow);

    void readTriangleModelData(const QString& filename, bool skipFirstRow);

    void readStandCatalog(const QString& filename, Catalog& catalog);

    void readRealData(const QString& filename, bool skipFirstRow, bool reverse = false);

    void readLines(const QString& filename, bool skipFirstRow, Catalog& catalog);

    void readLinesWithCatalog(const QString& filename, bool skipFirstRow, Catalog& catalog);

    void calculateFull(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors);

    void calculateTriangle(bool distorsio, Results& results, ResultErrors& errors);

    void calculateOldModel(const Catalog& catalog, OldErrors errors, double focus, qint32 numP);

    void fitFocusByLines(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors);

    void setPixelSize(double _pixelSize) {pixelSize = _pixelSize;}

    void setFrameSize(quint32 x, quint32 y) {frameX = x; frameY = y;}

    void setMeasureTheshold(quint32 th){thershold = th;}

    void findDistorsio(int nPow, bool clear = false);

    void saveDistorsio();

    void saveShifts(const QString& prefix);

    void includeDistorsio();

    QString printTestTable(const QString& filename, bool dist, double focus);

    QList <double> getDistX() const {return xDistV;}

    QList <double> getDistY() const {return yDistV;}

    ShiftData getShiftData() const;

    void setAxisDirection(X_AXIS_DIRECTION dir) {d = dir;}

    void clearAll(){frame.clear(); xDistV.clear(); yDistV.clear();
                    distData.x.clear(); distData.y.clear(); distData.dx.clear(); distData.dy.clear();}

    void setModelName(const QString& name) {modelName = name;}

    bool fishEye = false;

private:
    void includeAxisDirection (double MStand[3][3], double modifyMStand[3][3], X_AXIS_DIRECTION d = X_AXIS_DIRECTION::UP);

    void calculatePrivate(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors, const RotateAngles& rotAngles, QVector<QPointF>& frame);

    void findDistCft(int Npow, QVector <double>& x, QVector <double>& y, QVector <double>& dx, QVector <double>& dy);

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

    int gaussObr(int cntStar, double mass[55][55], double mObr[55][55]);

    double calculateAngle(const QPointF& fPoint, const QPointF& sPoint, double focus, QList <double>& distorsioCoefX, QList <double>& distorsioCoefY);

    void firstApprox(const Catalog& catalog, const LineData& lineData, OrientAngles& angles);

    void calculateOldModelPrivate(const Catalog& catalog, const LineData& lineData, OldErrors errors, double& focus, qint32 numP, bool dist, bool save);


    double pixelSize;
    quint32 thershold = 2;
    quint32 frameX;
    quint32 frameY;
    QVector <QPointF> frame;
    QVector <QPointF> initFrame;
    QVector < QVector <TrianglePattern> > triangleFrame;
    QVector <LineData> lineFrame;
    RotateAngles rotAngles;
    DistorsioData distData;
    QList <double> xDistV;
    QList <double> yDistV;
    X_AXIS_DIRECTION d = UP;
    QString modelName;

    static constexpr const double deltaAngle = 1. / 60. / 60. * degreesToRad * 10;
    static constexpr const double deltaFocus = 0.0001;
    static constexpr const qint32 maxDer = 10000;
    static constexpr const qint32 maxParams = 55;
    static constexpr const qint32 maxIterations = 500;
    static constexpr const double maxError = 0.000005;
};

#endif // MLSTASK_H
