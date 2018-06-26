#include "mlstask.h"
#include <QDebug>
MLSTask::MLSTask(QObject *parent) : QObject(parent)
{

}


void MLSTask::readModelData(const QString& filename, bool skipFirstRow)
{
    frame.clear();
    rotAngles.alphaRotates.clear();
    rotAngles.phiRotates.clear();
    QFile file (filename);
    if(file.open(QIODevice::ReadOnly))
    {
        QTextStream in(&file);
        QString line;
        constexpr const quint16 rowLength = 8;
        if (skipFirstRow)
        {
            in.readLineInto(&line);
        }
        while(in.readLineInto(&line))
        {
            QStringList list = line.split("\t");
            if (list.size() != rowLength)
                throw std::logic_error("Неверная строка в каталоге");
            frame.append(
                        QPointF((list[5].toDouble() - frameX) * pixelSize,
                        (list[6].toDouble() - frameY) * pixelSize));
            rotAngles.alphaRotates.append(list[3].toDouble());
            rotAngles.phiRotates.append(list[2].toDouble());
        }

        file.close();
    }
    else
    {
        throw std::runtime_error("Не удалось открыть файл " + filename.toStdString());
    }
}

void MLSTask::readRealData(const QString& filename, bool skipFirstRow)
{
    frame.clear();
    rotAngles.alphaRotates.clear();
    rotAngles.phiRotates.clear();
    QFile file (filename);
    if(file.open(QIODevice::ReadOnly))
    {
        QTextStream in(&file);
        QString line;
        constexpr const quint16 rowLength = 4;
        if (skipFirstRow)
        {
            in.readLineInto(&line);
        }
        while(in.readLineInto(&line))
        {
            QStringList list = line.split("\t");
            if (list.size() != rowLength)
                throw std::logic_error("Неверная строка в каталоге");
            frame.append(
                        QPointF((list[0].toDouble() - frameX) * pixelSize,
                        (list[1].toDouble() - frameY) * pixelSize));
            double alpha = list[2].toDouble();
            alpha = 95.0 + alpha;
            alpha = -alpha;
            double azimut = list[3].toDouble();
            //            if (frameX != frameY)
            //                azimut += 180.0;
            rotAngles.alphaRotates.append(alpha);
            rotAngles.phiRotates.append(azimut);
        }

        file.close();
    }
    else
    {
        throw std::runtime_error("Не удалось открыть файл " + filename.toStdString());
    }
}

void MLSTask::calculateXY(double Mstand[3][3], double collimator[3], double focus, double& X, double& Y)
{
    double LMN[3];
    LMN[0] = Mstand[0][0] * collimator[0] + Mstand[0][1] * collimator[1] + Mstand[0][2] * collimator[2];
    LMN[1] = Mstand[1][0] * collimator[0] + Mstand[1][1] * collimator[1] + Mstand[1][2] * collimator[2];
    LMN[2] = Mstand[2][0] * collimator[0] + Mstand[2][1] * collimator[1] + Mstand[2][2] * collimator[2];
    X = -focus * LMN[0] / LMN[2];
    Y = -focus * LMN[1] / LMN[2];
}


void MLSTask::includeAxisDirection (double MStand[3][3], double modifyMStand[3][3], X_AXIS_DIRECTION d)
{
    double m[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
    switch (d)
    {
    case UP:
    {
        m[0][0] = 1;
        m[1][1] = 1;
        m[2][2] = 1;
        break;
    }
    case DOWN:
    {
        m[0][0] = -1;
        m[1][1] = -1;
        m[2][2] = 1;
        break;
    }
    case LEFT:
    {
        m[0][1] = -1;
        m[1][0] = 1;
        m[2][2] = 1;
        break;
    }
    case RIGHT:
    {
        m[0][1] = 1;
        m[1][0] = -1;
        m[2][2] = 1;
        break;
    }

    }
    BOKZMath::multiplyMatrix(m, MStand, modifyMStand);
}

double MLSTask::calculate11El(StandAngles ang, double alpha, double phi)
{
    return cos(ang.gammaOY) * (sin(alpha) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) - cos(alpha) * sin(ang.alphaTwoOY) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            sin(ang.gammaOY) * (cos(ang.alphaTwoOY) * cos(alpha));
}

double MLSTask::calculate21El(StandAngles ang, double alpha, double phi)
{

    return sin(ang.gammaOY) * sin(ang.gammaOX) * (sin(alpha) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) - cos(alpha) * sin(ang.alphaTwoOY) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            cos(ang.gammaOX) * (-(cos(alpha) * sin(ang.alphaTwoOY) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi))) + sin(alpha) * (-(cos(ang.alphaTwoOX) * cos(phi)) + sin(ang.alphaTwoOX) * sin(phi)))
            +
            -(cos(ang.gammaOY) * sin(ang.gammaOX)) * (cos(ang.alphaTwoOY) * cos(alpha));
}

double MLSTask::calculate31El(StandAngles ang, double alpha, double phi)
{
    return -(cos(ang.gammaOX) * sin(ang.gammaOY)) * (sin(alpha) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) - cos(alpha) * sin(ang.alphaTwoOY) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            sin(ang.gammaOX) * (-(cos(alpha) * sin(ang.alphaTwoOY) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi))) + sin(alpha) * (-(cos(ang.alphaTwoOX) * cos(phi)) + sin(ang.alphaTwoOX) * sin(phi)))
            +
            cos(ang.gammaOY) * cos(ang.gammaOX) * (cos(ang.alphaTwoOY) * cos(alpha));
}

double MLSTask::calculate12El(StandAngles ang, double alpha, double phi)
{
    return cos(ang.gammaOY) * (cos(alpha) * cos(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) + cos(ang.lambdaOX) * sin(ang.alphaTwoOY) * sin(alpha) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)) + cos(ang.alphaTwoOY) * sin(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            sin(ang.gammaOY) * (-(cos(ang.alphaTwoOY) * cos(ang.lambdaOX) * sin(alpha)) + sin(ang.alphaTwoOY) * sin(ang.lambdaOX));
}

double MLSTask::calculate22El(StandAngles ang, double alpha, double phi)
{
    return sin(ang.gammaOY) * sin(ang.gammaOX) * (cos(alpha) * cos(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) + cos(ang.lambdaOX) * sin(ang.alphaTwoOY) * sin(alpha) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)) + cos(ang.alphaTwoOY) * sin(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            cos(ang.gammaOX) * (cos(ang.lambdaOX) * sin(ang.alphaTwoOY) * sin(alpha) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) + cos(ang.alphaTwoOY) * sin(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) + cos(alpha) * cos(ang.lambdaOX) * (-(cos(ang.alphaTwoOX) * cos(phi)) + sin(ang.alphaTwoOX) * sin(phi)))
            +
            -(cos(ang.gammaOY) * sin(ang.gammaOX)) * (-(cos(ang.alphaTwoOY) * cos(ang.lambdaOX) * sin(alpha)) + sin(ang.alphaTwoOY) * sin(ang.lambdaOX));
}

double MLSTask::calculate32El(StandAngles ang, double alpha, double phi)
{
    return -(cos(ang.gammaOX) * sin(ang.gammaOY)) * (cos(alpha) * cos(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) + cos(ang.lambdaOX) * sin(ang.alphaTwoOY) * sin(alpha) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)) + cos(ang.alphaTwoOY) * sin(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            sin(ang.gammaOX) * (cos(ang.lambdaOX) * sin(ang.alphaTwoOY) * sin(alpha) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) + cos(ang.alphaTwoOY) * sin(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) + cos(alpha) * cos(ang.lambdaOX) * (-(cos(ang.alphaTwoOX) * cos(phi)) + sin(ang.alphaTwoOX) * sin(phi)))
            +
            cos(ang.gammaOY) * cos(ang.gammaOX) * (-(cos(ang.alphaTwoOY) * cos(ang.lambdaOX) * sin(alpha)) + sin(ang.alphaTwoOY) * sin(ang.lambdaOX));
}

double MLSTask::calculate13El(StandAngles ang, double alpha, double phi)
{
    return cos(ang.gammaOY) * (-(cos(alpha) * sin(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi))) + cos(ang.alphaTwoOY) * cos(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)) - sin(ang.alphaTwoOY) * sin(alpha) * sin(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            sin(ang.gammaOY) * (cos(ang.lambdaOX) * sin(ang.alphaTwoOY) + cos(ang.alphaTwoOY) * sin(alpha) * sin(ang.lambdaOX));
}

double MLSTask::calculate23El(StandAngles ang, double alpha, double phi)
{
    return sin(ang.gammaOY) * sin(ang.gammaOX) * (-(cos(alpha) * sin(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi))) + cos(ang.alphaTwoOY) * cos(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)) - sin(ang.alphaTwoOY) * sin(alpha) * sin(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            cos(ang.gammaOX) * (cos(ang.alphaTwoOY) * cos(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) - sin(ang.alphaTwoOY) * sin(alpha) * sin(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) - cos(alpha) * sin(ang.lambdaOX) * (-(cos(ang.alphaTwoOX) * cos(phi)) + sin(ang.alphaTwoOX) * sin(phi)))
            +
            -(cos(ang.gammaOY) * sin(ang.gammaOX)) * (cos(ang.lambdaOX) * sin(ang.alphaTwoOY) + cos(ang.alphaTwoOY) * sin(alpha) * sin(ang.lambdaOX));

}

double MLSTask::calculate33El(StandAngles ang, double alpha, double phi)
{
    return -(cos(ang.gammaOX) * sin(ang.gammaOY)) * (-(cos(alpha) * sin(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi))) + cos(ang.alphaTwoOY) * cos(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)) - sin(ang.alphaTwoOY) * sin(alpha) * sin(ang.lambdaOX) * (cos(ang.alphaTwoOX) * cos(phi) - sin(ang.alphaTwoOX) * sin(phi)))
            +
            sin(ang.gammaOX) * (cos(ang.alphaTwoOY) * cos(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) - sin(ang.alphaTwoOY) * sin(alpha) * sin(ang.lambdaOX) * (cos(phi) * sin(ang.alphaTwoOX) + cos(ang.alphaTwoOX) * sin(phi)) - cos(alpha) * sin(ang.lambdaOX) * (-(cos(ang.alphaTwoOX) * cos(phi)) + sin(ang.alphaTwoOX) * sin(phi)))
            +
            cos(ang.gammaOY) * cos(ang.gammaOX) * (cos(ang.lambdaOX) * sin(ang.alphaTwoOY) + cos(ang.alphaTwoOY) * sin(alpha) * sin(ang.lambdaOX));
}

void MLSTask::calculateMatrix(StandAngles ang, double alpha, double phi, double Mstand[3][3])
{
    Mstand[0][0] = calculate11El(ang, alpha, phi);

    Mstand[1][0] = calculate21El(ang, alpha, phi);

    Mstand[2][0] = calculate31El(ang, alpha, phi);


    Mstand[0][1] = calculate12El(ang, alpha, phi);

    Mstand[1][1] = calculate22El(ang, alpha, phi);

    Mstand[2][1] = calculate32El(ang, alpha, phi) ;


    Mstand[0][2] = calculate13El(ang, alpha, phi);

    Mstand[1][2] = calculate23El(ang, alpha, phi);

    Mstand[2][2] = calculate33El(ang, alpha, phi);
}

void MLSTask::calculateMatrixDynamicaly(StandAngles ang, double alpha, double phi, double Mstand[3][3])
{
    double Minit[3][3];
    Minit[0][0] = 1; Minit[0][1] = 0; Minit[0][2] = 0;
    Minit[1][0] = 0; Minit[1][1] = 1; Minit[1][2] = 0;
    Minit[2][0] = 0; Minit[2][1] = 0; Minit[2][2] = 1;

    // матрица СК2-ВСК
    double MgammaOZ[3][3];
    rotateOZ(ang.gammaOZ, Minit, MgammaOZ);
    double MgammaOY[3][3];
    rotateOY(ang.gammaOY, MgammaOZ, MgammaOY);
    double MgammaOX[3][3];
    rotateOX(ang.gammaOX, MgammaOY, MgammaOX);

    // матрица ССК-СК1
   // double MfirstPlatfOX[3][3];
   //rotateOX(ang.lambdaOX, Minit, MfirstPlatfOX);
    double Mx[3][3];
    rotateOZ(alpha, Minit, Mx);

    // матрица СК2-ВСК
    double MsecondPlatfOY[3][3];
    rotateOY(ang.alphaTwoOY, Minit, MsecondPlatfOY);
    double MsecondPlatfOX[3][3];
    rotateOX(ang.alphaTwoOX, MsecondPlatfOY, MsecondPlatfOX);

    double Mtemp[3][3];
    double Mi[3][3]  = {{0, 0, 1}, {0, -1, 0}, {1, 0, 0}};
    multiplyMatrix(Mi, MsecondPlatfOX, Mtemp);
    double Mc [3][3];
    rotateOZ(phi, Mtemp, Mc);
    // итоговое перемножение
    multiplyMatrix(Mc ,Mx, Mtemp);
    multiplyMatrix(MgammaOX, Mtemp, Mstand);

}

void MLSTask::fitFocusByLines(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors)
{

    QVector <double> focs;
    bool lastLine = false;
    RotateAngles angs;
    QVector<QPointF> smallFrame;
    int curPos = 0;
    int size = rotAngles.phiRotates.size();
    while(!lastLine)
    {
        double angle = rotAngles.phiRotates[curPos];
        while (curPos < size
               && qFuzzyCompare(rotAngles.phiRotates[curPos], angle))
        {
            angs.alphaRotates.append(rotAngles.alphaRotates[curPos]);
            angs.phiRotates.append(rotAngles.phiRotates[curPos]);
            smallFrame.append(QPointF(frame[curPos].x(),frame[curPos].y()));
            curPos++;
        }
        if (curPos == size)
        {
            results.foc = calculateMean(focs.begin(), focs.end(), 0.0);
            lastLine = true;
            continue;
        }
        calculatePrivate(derivativeFlags,results,errors, angs, smallFrame);
        focs.append(results.foc);
        angs.alphaRotates.clear();
        angs.phiRotates.clear();
        smallFrame.clear();
    }
}
void MLSTask::calculate(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors)
{
    calculatePrivate(derivativeFlags, results, errors, rotAngles, frame);

}
void MLSTask::calculatePrivate(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors, const RotateAngles& rotAngles, QVector<QPointF>& frame)
{
    double MstandTemp[3][3];
    double Mstand[3][3];
    double MstandDer[3][3];
    StandAngles ang;
    StandAngles angDer;

    ang.gammaOZ = results.gammaOZ * degreesToRad;
    ang.gammaOY = results.gammaOY * degreesToRad;
    ang.gammaOX = results.gammaOX * degreesToRad;
    ang.lambdaOY = results.lambdaOY * degreesToRad;
    ang.lambdaOX = results.lambdaOX * degreesToRad;
    ang.alphaTwoOY = results.alphaTwoOY * degreesToRad;
    ang.alphaTwoOX = results.alphaTwoOX * degreesToRad;


    double collimator[3] {-(cos(ang.lambdaOY) * cos(ang.lambdaOX)), -(cos(ang.lambdaOY) * sin(ang.lambdaOX)), -sin(ang.lambdaOY)};
    double collimatorDer[3];

    QVector <double> Xst;
    QVector <double> Yst;
    for (int i = 0; i < frame.size(); i++)
    {
        Xst.append(frame[i].x());
        Yst.append(frame[i].y());
    }

    double DFR[maxDer];
    double **DRV = new double* [maxDer];
    double **DRVM = new double* [maxParams];
    double **DRVM_1 = new double* [maxParams];
    for (int count = 0; count < maxDer; count++)
    {
        DRV[count] = new double [maxParams];
    }
    for (int j = 0; j < maxParams; j++)
    {
        DRVM[j] = new double [maxParams];
        DRVM_1[j] = new double [maxParams];
    }

    double DRVH[maxParams], ANGD[maxParams];
    double mx, my;
    double mxy = 1000;
    double mxy_pr;
    qint32 clk = 0;
    qint32 Nst = frame.size();
    double X, Y, Xd, Yd;
    qint32 numP = 0;
    double focus = results.foc;
    do
    {
        mxy_pr = mxy;
        mxy = 0.;  mx = 0.; my = 0.;

        for (int j = 0; j < Nst; j++)
        {
            numP = 0;
            int k = j * 2; // номер уравнения для х
            int l = k + 1; // номер уравнения для у

            double alpha = rotAngles.alphaRotates[j] * degreesToRad;
            double phi = rotAngles.phiRotates[j] * degreesToRad;
            calculateMatrixDynamicaly(ang, alpha, phi, MstandTemp);
            includeAxisDirection(MstandTemp, Mstand, d);
            calculateXY(Mstand, collimator, focus, X, Y);

            DFR[k] = (double)(Xst[j] - X);
            DFR[l] = (double)(Yst[j] - Y);


            if (derivativeFlags.at(DERIVATIVES::LAMBDA_OY))
            {
                angDer = ang;
                angDer.lambdaOY += deltaAngle;
                collimatorDer[0] = -(cos(angDer.lambdaOY) * cos(angDer.lambdaOX));
                collimatorDer[1] = -(cos(angDer.lambdaOY) * sin(angDer.lambdaOX));
                collimatorDer[2] = -sin(angDer.lambdaOY);
                calculateMatrixDynamicaly(angDer, alpha, phi, MstandTemp);
                includeAxisDirection(MstandTemp, MstandDer, d);
                calculateXY(MstandDer, collimatorDer, focus, Xd, Yd);
                DRV[k][numP] = (Xd - X) / deltaAngle;
                DRV[l][numP++] = (Yd - Y) / deltaAngle;
            }


            if (derivativeFlags.at(DERIVATIVES::LAMBDA_OX))
            {
                angDer = ang;
                angDer.lambdaOX += deltaAngle;
                collimatorDer[0] = -(cos(angDer.lambdaOY) * cos(angDer.lambdaOX));
                collimatorDer[1] = -(cos(angDer.lambdaOY) * sin(angDer.lambdaOX));
                collimatorDer[2] = -sin(angDer.lambdaOY);
                calculateMatrixDynamicaly(angDer, alpha, phi, MstandTemp);
                includeAxisDirection(MstandTemp, MstandDer, d);
                calculateXY(MstandDer, collimatorDer, focus, Xd, Yd);
                DRV[k][numP] = (Xd - X) / deltaAngle;
                DRV[l][numP++] = (Yd - Y) / deltaAngle;
            }

            if (derivativeFlags.at(DERIVATIVES::ALPHA_TWO_OY))
            {
                angDer = ang;
                angDer.alphaTwoOY += deltaAngle;
                calculateMatrixDynamicaly(angDer, alpha, phi, MstandTemp);
                includeAxisDirection(MstandTemp, MstandDer, d);
                calculateXY(MstandDer, collimator, focus, Xd, Yd);
                DRV[k][numP] = (Xd - X) / deltaAngle;
                DRV[l][numP++] = (Yd - Y) / deltaAngle;
            }

            if (derivativeFlags.at(DERIVATIVES::ALPHA_TWO_OX))
            {
                angDer = ang;
                angDer.alphaTwoOX += deltaAngle;
                calculateMatrixDynamicaly(angDer, alpha, phi, MstandTemp);
                includeAxisDirection(MstandTemp, MstandDer, d);
                calculateXY(MstandDer, collimator, focus, Xd, Yd);
                DRV[k][numP] = (Xd - X) / deltaAngle;
                DRV[l][numP++] = (Yd - Y) / deltaAngle;
            }

            if (derivativeFlags.at(DERIVATIVES::GAMMA_OZ))
            {
                angDer = ang;
                angDer.gammaOZ += deltaAngle;
                calculateMatrixDynamicaly(angDer, alpha, phi, MstandTemp);
                includeAxisDirection(MstandTemp, MstandDer, d);
                calculateXY(MstandDer, collimator, focus, Xd, Yd);
                DRV[k][numP] = (Xd - X) / deltaAngle;
                DRV[l][numP++] = (Yd - Y) / deltaAngle;
            }

            if (derivativeFlags.at(DERIVATIVES::GAMMA_OY))
            {
                angDer = ang;
                angDer.gammaOY += deltaAngle;
                calculateMatrixDynamicaly(angDer, alpha, phi, MstandTemp);
                includeAxisDirection(MstandTemp, MstandDer, d);
                calculateXY(MstandDer, collimator, focus, Xd, Yd);
                DRV[k][numP] = (Xd - X) / deltaAngle;
                DRV[l][numP++] = (Yd - Y) / deltaAngle;
            }


            if (derivativeFlags.at(DERIVATIVES::GAMMA_OX))
            {
                angDer = ang;
                angDer.gammaOX += deltaAngle;
                calculateMatrixDynamicaly(angDer, alpha, phi, MstandTemp);
                includeAxisDirection(MstandTemp, MstandDer, d);
                calculateXY(MstandDer, collimator, focus, Xd, Yd);
                DRV[k][numP] = (Xd - X) / deltaAngle;
                DRV[l][numP++] = (Yd - Y) / deltaAngle;

            }
            if (derivativeFlags.at(DERIVATIVES::FOCUS))
            {
                calculateXY(Mstand, collimator, focus + deltaFocus, Xd, Yd);
                DRV[k][numP] = (Xd - X) / deltaFocus;
                DRV[l][numP++] = (Yd - Y) / deltaFocus;
            }

            mxy += DFR[k] * DFR[k] + DFR[l] * DFR[l];
            mx  += DFR[k] * DFR[k];
            my  += DFR[l] * DFR[l];
        }
        mxy = mxy / (double)(2 * Nst - numP);  mxy = sqrtm(mxy);
        mx  = mx / (double)(Nst - numP/2.);    mx  = sqrtm(mx);
        my  = my / (double)(Nst - numP/2.);    my  = sqrtm(my);

        for (int i = 0; i < numP; i++)
        {
            for (int j = 0; j < numP; j++) // произведение обычной и транс. матрицы параметров
            {
                DRVM[i][j] = 0.;      //(DRV)^T*DRV
                for (int k = 0; k < 2 * Nst; k++)
                {
                    DRVM[i][j] += DRV[k][i] * DRV[k][j];
                }
            }
        }
        if (numP == 0)
        {
            return;
        }
        else if (numP == 1)
        {
            DRVM_1[0][0] = 1 / DRVM[0][0];
        }
        else
        {
            Matrix_1MM(DRVM, DRVM_1, numP); // вычисление обратной матрицы
        }
        for (int i = 0; i < numP; i++)
        { // умножение транспонированной на вектор расхождений координат
            DRVH[i] = 0.;             //DRVH=(DRV)^T*DIF
            for (int k = 0; k < 2 * Nst; k++)
            {
                DRVH[i] += DRV[k][i] * DFR[k];
            }
        }

        for (int i = 0; i < numP; i++)
        { // вычисление дельт углов
            ANGD[i] = 0.;
            for (int k = 0; k < numP; k++)
            {
                ANGD[i] += DRVM_1[i][k] * DRVH[k];
            }
        }

        numP = 0;
        if (derivativeFlags.at(DERIVATIVES::LAMBDA_OY))
        {
            ang.lambdaOY += ANGD[numP++];
        }

        if (derivativeFlags.at(DERIVATIVES::LAMBDA_OX))
        {
            ang.lambdaOX += ANGD[numP++];
        }

        if (derivativeFlags.at(DERIVATIVES::ALPHA_TWO_OY))
        {
            ang.alphaTwoOY += ANGD[numP++];
        }

        if (derivativeFlags.at(DERIVATIVES::ALPHA_TWO_OX))
        {
            ang.alphaTwoOX += ANGD[numP++];
        }

        if (derivativeFlags.at(DERIVATIVES::GAMMA_OZ))
        {
            ang.gammaOZ += ANGD[numP++];
        }

        if (derivativeFlags.at(DERIVATIVES::GAMMA_OY))
        {
            ang.gammaOY += ANGD[numP++];
        }

        if (derivativeFlags.at(DERIVATIVES::GAMMA_OX))
        {
            ang.gammaOX += ANGD[numP++];
        }

        if (derivativeFlags.at(DERIVATIVES::FOCUS))
        {
            focus += ANGD[numP++];
        }

//        if(fabs((ang.lambdaOY * radToDegrees) - results.lambdaOY) > 2)
//            ang.lambdaOY = results.lambdaOY * degreesToRad;

//        if(fabs((ang.lambdaOX * radToDegrees) - results.lambdaOX)  > 10)
//            ang.lambdaOX = results.lambdaOX * degreesToRad;

//        if(fabs((ang.alphaTwoOY * radToDegrees) - results.alphaTwoOY) > 2)
//            ang.alphaTwoOY = results.alphaTwoOY * degreesToRad;

//        if(fabs((ang.alphaTwoOX * radToDegrees) - results.alphaTwoOX) > 2)
//            ang.alphaTwoOX = results.alphaTwoOX * degreesToRad;

//        if(fabs((ang.gammaOZ * radToDegrees) - results.gammaOZ) > 2)
//            ang.gammaOZ = results.gammaOZ * degreesToRad;

//        if(fabs((ang.gammaOY * radToDegrees) - results.gammaOY) > 6)
//            ang.gammaOY = results.gammaOY * degreesToRad;

//        if(fabs((ang.gammaOX * radToDegrees)- results.gammaOX) > 6)
//            ang.gammaOX = results.gammaOX *  degreesToRad;

//        if(fabs(focus - results.foc) > 1) focus = results.foc;

        //это тоже не влияет

        collimator[0] = -(cos(ang.lambdaOY) * cos(ang.lambdaOX));
        collimator[1] = -(cos(ang.lambdaOY) * sin(ang.lambdaOX));
        collimator[2] = -sin(ang.lambdaOY);

//        collimator[0] = -cos(ang.lambdaOY);
//        collimator[1] = 0;
//        collimator[2] = -sin(ang.lambdaOY);
        clk++;
    }
    while ((fabs(mxy - mxy_pr) > maxError)
           && (clk < maxIterations));



    results.foc = focus;
    results.alphaTwoOX = ang.alphaTwoOX * radToDegrees;
    results.alphaTwoOY = ang.alphaTwoOY * radToDegrees;
    results.lambdaOY = ang.lambdaOY * radToDegrees;
    results.lambdaOX = ang.lambdaOX * radToDegrees;
    results.gammaOZ = ang.gammaOZ * radToDegrees;
    results.gammaOY = ang.gammaOY * radToDegrees;
    results.gammaOX = ang.gammaOX * radToDegrees;

    numP = 0;
    if (derivativeFlags.at(DERIVATIVES::LAMBDA_OY))
    {
        errors.dlambdaOY = sqrtm(DRVM_1[numP][numP]) * mxy * radToSec;
        ++numP;
    }

    if (derivativeFlags.at(DERIVATIVES::LAMBDA_OX))
    {
        errors.dlambdaOX = sqrtm(DRVM_1[numP][numP]) * mxy * radToSec;
        ++numP;
    }

    if (derivativeFlags.at(DERIVATIVES::ALPHA_TWO_OY))
    {
        errors.dalphaTwoOY = sqrtm(DRVM_1[numP][numP]) * mxy * radToSec;
        ++numP;
    }

    if (derivativeFlags.at(DERIVATIVES::ALPHA_TWO_OX))
    {
        errors.dalphaTwoOX = sqrtm(DRVM_1[numP][numP]) * mxy * radToSec;
        ++numP;
    }

    if (derivativeFlags.at(DERIVATIVES::GAMMA_OZ))
    {
        errors.dgammaOZ = sqrtm(DRVM_1[numP][numP]) * mxy * radToSec;
        ++numP;
    }

    if (derivativeFlags.at(DERIVATIVES::GAMMA_OY))
    {
        errors.dgammaOY = sqrtm(DRVM_1[numP][numP]) * mxy * radToSec;
        ++numP;
    }

    if (derivativeFlags.at(DERIVATIVES::GAMMA_OX))
    {
        errors.dgammaOX = sqrtm(DRVM_1[numP][numP]) * mxy * radToSec;
        ++numP;
    }
    if (derivativeFlags.at(DERIVATIVES::FOCUS))
    {
        errors.dfoc = sqrtm(DRVM_1[numP][numP]) * mxy * 1000;
        ++numP;
    }

    errors.mx = mx * 1000;
    errors.my = my * 1000;
    errors.mxy = mxy * 1000;


    if (mxy > thershold * pixelSize)
    {
        throw std::runtime_error("Грубые измерения.");
    }

    distData.x.clear();
    distData.y.clear();
    distData.dx.clear();
    distData.dy.clear();
    //qDebug() << "рассогласования";
    for (int i = 0; i < Nst; i++)
    {
        int deltaStep = i * 2;
        distData.x.append(Xst[i]);
        distData.y.append(Yst[i]);
        distData.dx.append(DFR[deltaStep]);
        distData.dy.append(DFR[deltaStep + 1]);
    }

    for (int count = 0; count < Nst; count++)
    {
        delete[] DRV[count];
    }
    delete[] DRV;
    for (int count = 0; count < maxParams; count++)
    {
        delete[] DRVM[count];
        delete[] DRVM_1[count];
    }
    delete [] DRVM;
    delete [] DRVM_1;
}


void MLSTask::findDistorsio(int nPow)
{
    FindDistCft(nPow, distData.x, distData.y, distData.dx, distData.dy);
}


void MLSTask::saveDistorsio()
{
    QFile file("dist.txt");
    if (file.open(QIODevice::WriteOnly))
    {
        QTextStream out(&file);
        for (int i = 0; i < xDistV.size(); i++)
        {
            out << QString("%1 %2").arg(xDistV[i], 0, 'f', 15).arg(yDistV[i], 0, 'f', 15) << "\n";
        }
    }
}

void MLSTask::saveShifts(const QString& prefix)
{
    QFile coords(prefix + "coodinates" + ".txt");
    if (coords.open(QIODevice::WriteOnly))
    {
        QTextStream out(&coords);
        for (int i = 0; i < distData.x.size(); i++)
        {
            out << QString("%1 %2 %3 %4 %5 %6 %7 %8\n")
                   .arg(frame[i].x())
                   .arg(frame[i].y())
                   .arg(distData.x[i])
                   .arg(distData.y[i])
                   .arg(distData.dx[i])
                   .arg(distData.dy[i])
                   .arg(distData.dx_diff[i])
                   .arg(distData.dy_diff[i]);
        }
    }
}

void MLSTask::includeDistorsio()
{

    for (int i = 0; i < frame.size(); i++)
    {
        frame[i].setX(calculateDistorsio(frame[i].x(), frame[i].x(), frame[i].y(), xDistV));
        frame[i].setY(calculateDistorsio(frame[i].y(), frame[i].x(), frame[i].y(), yDistV));
    }
}

QVector <QString> MLSTask::printTestTable(const QString& filename, bool dist, double focus)
{
    double frameXF = frameX * 2;
    double frameYF = frameY * 2;
    QVector <QPointF> points {
        QPointF(0,0),
                QPointF(frameXF,0),
                QPointF(frameXF, frameYF),
                QPointF(0, frameYF),
                QPointF(frameXF * (1./2),frameYF * (1./3)),
                QPointF(frameXF * (2./3),frameYF * (1./2)),
                QPointF(frameXF * (1./2),frameYF * (2./3)),
                QPointF(frameXF * (1./3),frameYF * (1./2))};

    for(int i = 0; i < points.size(); i++)
    {
        points[i].setX((points[i].x() - frameX) * pixelSize);
        points[i].setY((points[i].y() - frameY) * pixelSize);
    }
    if (dist)
    {
        for (int i = 0; i < points.size(); i++)
        {
            points[i].setX(calculateDistorsio(points[i].x(), points[i].x(), points[i].y(), xDistV));
            points[i].setY(calculateDistorsio(points[i].y(), points[i].x(), points[i].y(), yDistV));
        }
    }

    QVector <QVector <double>> lmn;
    for (int i = 0; i < points.size(); i++)
    {
        double lmnt[3];
        BOKZMath::calculateLMNImage(points[i].x(), points[i].y(), focus, lmnt);
        lmn.append(QVector <double> {lmnt[0], lmnt[1], lmnt[2]});
    }
    double m[8][8];
    for (int i = 0; i < points.size(); i++)
    {
        for (int j = 0; j < points.size(); j++)
        {
            double ac = acos(calculateScalarProduct
                             (lmn[i][0], lmn[j][0], lmn[i][1], lmn[j][1], lmn[i][2], lmn[j][2]));
            m[i][j] = ac * radToDegrees;
        }
    }

    QString line;
    QString retLine;
    QVector <QString> retLineVec;
    for (int i = 0; i < points.size(); i++)
    {
        line.clear();
        for (int j = 0; j < points.size(); j++)
        {
            GMS g(m[i][j]);
            line.append( QString("%1 %2'%3''      ")
                         .arg(g.getGradus())
                         .arg(g.getMinutes())
                         .arg(g.getSeconds()));
            retLine.append(QString::number(m[i][j], 'f', 7) + "  ");

        }
        retLineVec.append(retLine);
        retLine.clear();
    }

    QFile f(filename);
    if (f.open(QIODevice::WriteOnly))
    {
        QTextStream out(&f);
        for (int i = 0; i < points.size(); i++)
        {
            line.clear();
            for (int j = 0; j < points.size(); j++)
            {
                GMS g(m[i][j]);
                out << QString("%1 %2'%3''\t")
                       .arg(g.getGradus())
                       .arg(g.getMinutes())
                       .arg(g.getSeconds());
            }
            out << "\n";
        }
    }
    return retLineVec;
}

void MLSTask::FindDistCft(int Npow, QVector <double>& x, QVector <double>& y, QVector <double>& dx, QVector <double>& dy)
{
    double CX[maxParams], CY[maxParams],B[maxParams][maxParams],
            B_1[maxParams][maxParams], AX[maxParams], AY[maxParams];
    double Mcft_i[maxParams];
    double x_2, y_2, x_3, y_3, x_4, y_4, x_5,
            y_5, x_6, y_6, x_7, y_7, x_8, y_8, x_9, y_9;
    int Ncft;

    //    QFile coords("coodinates" + number + ".txt");
    //    if (coords.open(QIODevice::WriteOnly))
    //    {
    //        QTextStream out(&coords);
    //        for (int i = 0; i < x.size(); i++)
    //        {
    //            out << QString("%1  %2  %3  %4\n").arg(x[i]).arg(y[i]).arg(dx[i]).arg(dy[i]);
    //        }
    //    }

    switch (Npow)
    {
    case 2: Ncft=6;  break;
    case 3: Ncft=10; break;
    case 4: Ncft=15; break;
    case 5: Ncft=21; break;
    case 6: Ncft=28; break;
    case 7: Ncft=36; break;
    case 8: Ncft=45; break;
    case 9: Ncft=55; break;
    default: Npow=3; Ncft=10; break;
    }


    for (int i = 0;i < maxParams; i++)
    {
        CX[i]=0.; CY[i]=0.;
        for (int j=0;j<Ncft;j++)
            B[i][j]=0.;
    }

    for (int i = 0; i < x.size(); i++)
    {
        x_2 = x[i] * x[i];
        y_2 = y[i] * y[i];

        Mcft_i[0] = 1.;
        Mcft_i[1] = x[i];
        Mcft_i[2] = y[i];
        Mcft_i[3] = x_2;
        Mcft_i[4] = x[i]*y[i];
        Mcft_i[5] = y_2;
        if (Npow>=3)
        {
            x_3=x_2*x[i];
            y_3=y_2*y[i];

            Mcft_i[6]=x_3;
            Mcft_i[7]=x_2*y[i];
            Mcft_i[8]=x[i]*y_2;
            Mcft_i[9]=y_3;
            if (Npow>=4)
            {
                x_4=x_3*x[i];
                y_4=y_3*y[i];

                Mcft_i[10]=x_4;
                Mcft_i[11]=x_3*y[i];
                Mcft_i[12]=x_2*y_2;
                Mcft_i[13]=x[i]*y_3;
                Mcft_i[14]=y_4;

                if (Npow>=5)
                {
                    x_5=x_4*x[i];
                    y_5=y_4*y[i];

                    Mcft_i[15]=x_5;
                    Mcft_i[16]=x_4*y[i];
                    Mcft_i[17]=x_3*y_2;
                    Mcft_i[18]=x_2*y_3;
                    Mcft_i[19]=x[i]*y_4;
                    Mcft_i[20]=y_5;
                    if (Npow>=6)
                    {
                        x_6=x_5*x[i];
                        y_6=y_5*y[i];

                        Mcft_i[21]=x_6;
                        Mcft_i[22]=x_5*y[i];
                        Mcft_i[23]=x_4*y_2;
                        Mcft_i[24]=x_3*y_3;
                        Mcft_i[25]=x_2*y_4;
                        Mcft_i[26]=x[i]*y_5;
                        Mcft_i[27]=y_6;
                        if (Npow>=7)
                        {
                            x_7=x_6*x[i];
                            y_7=y_6*y[i];

                            Mcft_i[28]=x_7;
                            Mcft_i[29]=x_6*y[i];
                            Mcft_i[30]=x_5*y_2;
                            Mcft_i[31]=x_4*y_3;
                            Mcft_i[32]=x_3*y_4;
                            Mcft_i[33]=x_2*y_5;
                            Mcft_i[34]=x[i]*y_6;
                            Mcft_i[35]=y_7;
                            if (Npow>=8)
                            {
                                x_8=x_7 * x[i];
                                y_8=y_7 * y[i];

                                Mcft_i[36]=x_8;
                                Mcft_i[37]=x_7 * y[i];
                                Mcft_i[38]=x_6*y_2;
                                Mcft_i[39]=x_5*y_3;
                                Mcft_i[40]=x_4*y_4;
                                Mcft_i[41]=x_3*y_5;
                                Mcft_i[42]=x_2*y_6;
                                Mcft_i[43]=x[i]*y_7;
                                Mcft_i[44]=y_8;
                                if (Npow >= 9)
                                {
                                    x_9=x_8*x[i];
                                    y_9=y_8*y[i];

                                    Mcft_i[45]=x_9;
                                    Mcft_i[46]=x_8*y[i];
                                    Mcft_i[47]=x_7*y_2;
                                    Mcft_i[48]=x_6*y_3;
                                    Mcft_i[49]=x_5*y_4;
                                    Mcft_i[50]=x_4*y_5;
                                    Mcft_i[51]=x_3*y_6;
                                    Mcft_i[52]=x_2*y_7;
                                    Mcft_i[53]=x[i]*y_8;
                                    Mcft_i[54]=y_9;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (int k = 0;k < Ncft; k++)
        {
            CX[k] += Mcft_i[k] * dx[i];
            CY[k] += Mcft_i[k] * dy[i];
            for (int j = 0;j < Ncft; j++)
                B[k][j] += Mcft_i[k] * Mcft_i[j];
        }

    }

    gaus_obr(Ncft,B,B_1);    //(A^T*A)^(-1)

    for (int i = Ncft; i < maxParams; i++)
    {
        AX[i] = 0.; AY[i] = 0;
    }

    for (int i = 0;i < Ncft; i++)
    {
        AX[i]=0.; AY[i]=0.;
        for (int j = 0;j < Ncft; j++)
        {
            AX[i] += B_1[i][j] * CX[j];     // (A^T*A)^(-1)*A^T*dx
            AY[i] += B_1[i][j] * CY[j];     // (A^T*A)^(-1)*A^T*dy
        }
    }
    xDistV.clear();
    yDistV.clear();
    constexpr const int maxDistCoef = 21;
    for (int i = 0; i < maxDistCoef; i++)
    {
        xDistV.append(0);
        yDistV.append(0);
    }
    for (int i = 0; i < Ncft; i++)
    {
        xDistV[i] = AX[i];
        yDistV[i] = AY[i];
    }


    //    QVector <double> mean_x;
    //    QVector <double> mean_y;
    //    QFile file("dist" + number + ".txt");
    //    if (file.open(QIODevice::WriteOnly))
    //    {
    //        QTextStream out(&file);
    //        for (int i = 0; i < Ncft; i++)
    //        {
    //            out << QString("%1 %2").arg(AX[i], 0, 'f', 15).arg(AY[i], 0, 'f', 15) << "\n";
    //        }
    //    }

    distData.dx_diff.clear();
    distData.dy_diff.clear();
    for (int i = 0; i < x.size(); i++)
    {

        x_2=x[i]*x[i];
        y_2=y[i]*y[i];
        x_3=x_2*x[i];
        y_3=y_2*y[i];
        x_4=x_3*x[i];
        y_4=y_3*y[i];
        x_5=x_4*x[i];
        y_5=y_4*y[i];
        x_6=x_5*x[i];
        y_6=y_5*y[i];
        x_7=x_6*x[i];
        y_7=y_6*y[i];
        x_8=x_7*x[i];
        y_8=y_7*y[i];
        x_9=x_8*x[i];
        y_9=y_8*y[i];

        double XP=AX[0]+AX[1]*x[i]+AX[2]*y[i]+AX[3]*x_2+AX[4]*x[i] * y[i]+AX[5]*y_2
                +AX[6]*x_3+AX[7]*x_2*y[i]+AX[8]*x[i]*y_2+AX[9]*y_3
                +AX[10]*x_4+AX[11]*x_3*y[i]+AX[12]*x_2*y_2+AX[13]*x[i]*y_3+AX[14]*y_4
                +AX[15]*x_5+AX[16]*x_4*y[i]+AX[17]*x_3*y_2+AX[18]*x_2*y_3+AX[19]*x[i]*y_4+AX[20]*y_5
                +AX[21]*x_6+AX[22]*x_5*y[i]+AX[23]*x_4*y_2+AX[24]*x_3*y_3+AX[25]*x_2*y_4+AX[26]*x[i]*y_5+AX[27]*y_6
                +AX[28]*x_7+AX[29]*x_6*y[i]+AX[30]*x_5*y_2+AX[31]*x_4*y_3+AX[32]*x_3*y_4+AX[33]*x_2*y_5+AX[34]*x[i]*y_6+AX[35]*y_7
                +AX[36]*x_8+AX[37]*x_7*y[i]+AX[38]*x_6*y_2+AX[39]*x_5*y_3+AX[40]*x_4*y_4+AX[41]*x_3*y_5+AX[42]*x_2*y_6+AX[43]*x[i]*y_7+AX[44]*y_8
                +AX[45]*x_9+AX[46]*x_8*y[i]+AX[47]*x_7*y_2+AX[48]*x_6*y_3+AX[49]*x_5*y_4+AX[50]*x_4*y_5+AX[51]*x_3*y_6+AX[52]*x_2*y_7+AX[53]*x[i]*y_8+AX[54]*y_9;

        double YP=AY[0]+AY[1]*x[i]+AY[2]*y[i]+AY[3]*x_2+AY[4]*x[i]*y[i]+AY[5]*y_2
                +AY[6]*x_3+AY[7]*x_2*y[i]+AY[8]*x[i]*y_2+AY[9]*y_3
                +AY[10]*x_4+AY[11]*x_3*y[i]+AY[12]*x_2*y_2+AY[13]*x[i]*y_3+AY[14]*y_4
                +AY[15]*x_5+AY[16]*x_4*y[i]+AY[17]*x_3*y_2+AY[18]*x_2*y_3+AY[19]*x[i]*y_4+AY[20]*y_5
                +AY[21]*x_6+AY[22]*x_5*y[i]+AY[23]*x_4*y_2+AY[24]*x_3*y_3+AY[25]*x_2*y_4+AY[26]*x[i]*y_5+AY[27]*y_6
                +AY[28]*x_7+AY[29]*x_6*y[i]+AY[30]*x_5*y_2+AY[31]*x_4*y_3+AY[32]*x_3*y_4+AY[33]*x_2*y_5+AY[34]*x[i]*y_6+AY[35]*y_7
                +AY[36]*x_8+AY[37]*x_7*y[i]+AY[38]*x_6*y_2+AY[39]*x_5*y_3+AY[40]*x_4*y_4+AY[41]*x_3*y_5+AY[42]*x_2*y_6+AY[43]*x[i]*y_7+AY[44]*y_8
                +AY[45]*x_9+AY[46]*x_8*y[i]+AY[47]*x_7*y_2+AY[48]*x_6*y_3+AY[49]*x_5*y_4+AY[50]*x_4*y_5+AY[51]*x_3*y_6+AY[52]*x_2*y_7+AY[53]*x[i]*y_8+AY[54]*y_9;


        double dx_new = dx[i] - XP;
        double dy_new = dy[i] - YP;
        distData.dx_diff.append(dx_new);
        distData.dy_diff.append(dy_new);
        //        qDebug() << dx_new << dy_new;
        // mean_x.append(fabs(dx_new));
        // mean_y.append(fabs(dy_new));

    }
}


int MLSTask::gaus_obr(int cnt_str,double mass[55][55],double M_obr[55][55])
{
    int i,j,k;
    double a,b;
    double sum;

    for(i=0;i<cnt_str;i++)
    {
        for(j=0;j<cnt_str;j++)
            M_obr[i][j]=0;
        M_obr[i][i]=1;
    }

    for(i=0;i<cnt_str;i++)
    {
        a=mass[i][i];
        if (fabs(a)<1e-38) return 0;
        for(j=i+1;j<cnt_str;j++)
        {
            b=mass[j][i];
            for(k=0;k<cnt_str;k++)
            {
                mass[j][k]=mass[j][k]-mass[i][k]*b/a;
                M_obr[j][k]=M_obr[j][k]-M_obr[i][k]*b/a;
            }
        }
    }

    for(i=0;i<cnt_str;i++)
    {
        for(j=cnt_str-1;j>=0;j--)
        {
            sum=0;
            for(k=cnt_str-1;k>j;k--)
                sum+=mass[j][k]*M_obr[k][i];
            if(mass[j][j]==0)
            {
                return 0;
            }
            M_obr[j][i]=(M_obr[j][i]-sum)/mass[j][j];
        }
    }
    return 1;
}
