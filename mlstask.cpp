#include "mlstask.h"
#include <QDebug>
MLSTask::MLSTask(QObject *parent) : QObject(parent)
{

}


void MLSTask::readFullModelData(const QString& filename, bool skipFirstRow)
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
        initFrame = frame;
        file.close();
    }
    else
    {
        throw std::runtime_error("Не удалось открыть файл " + filename.toStdString());
    }
}

void MLSTask::readTriangleModelData(const QString& filename, bool skipFirstRow)
{
    triangleFrame.clear();
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
        int pointsCount = 0;
        QString previousAlpha;
        QString previousAzimut;
        TrianglePattern triangle;
        in.readLineInto(&line);
        QStringList splitted = line.split("\t");

        previousAlpha = splitted[3];
        previousAzimut = splitted[2];
        initFrame.append(QPointF((splitted[5].toDouble() - frameX) * pixelSize,
                         (splitted[6].toDouble() - frameX) * pixelSize));
        triangle.points[pointsCount++] = initFrame.last();
        QVector <TrianglePattern> lineVec;
        while (in.readLineInto(&line))
        {
            QStringList list = line.split("\t");
            if (list.size() != rowLength)
            {
                throw std::logic_error("Неверная строка в каталоге");
            }
            if (previousAzimut != list[2])
            {
                triangleFrame.append(lineVec);
                previousAzimut = list[2];
                lineVec.clear();
            }
            if (previousAlpha != list[3])
            {
                previousAlpha = list[3];
                pointsCount = 0;
            }
            initFrame.append(QPointF((list[5].toDouble() - frameX) * pixelSize,
                             (list[6].toDouble() - frameY) * pixelSize));
            triangle.points[pointsCount++] = initFrame.last();

            if (pointsCount == countPointsTriangle)
            {
                lineVec.append(triangle);
                pointsCount = 0;
            }
        }
        file.close();
    }
    else
    {
        throw std::runtime_error("Не удалось открыть файл " + filename.toStdString());
    }
}

void MLSTask::readRealData(const QString& filename, bool skipFirstRow, bool reverse)
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
        while (in.readLineInto(&line))
        {
            QStringList list = line.split("\t");
            if (list.size() != rowLength)
                throw std::logic_error("Неверная строка в каталоге");

            if (!reverse)
            {
                frame.append(
                            QPointF((list[0].toDouble() - frameX) * pixelSize,
                            (list[1].toDouble() - frameY) * pixelSize));
            }
            else
            {
                frame.append(
                            QPointF((list[1].toDouble() - frameY) * pixelSize,
                            (list[0].toDouble() - frameX) * pixelSize));
            }

            double alpha = list[2].toDouble();
            alpha = 95.0 + alpha;
            alpha = -alpha;
            double azimut = list[3].toDouble();
            //            if (frameX != frameY)
            //azimut -= 180.0;
            rotAngles.alphaRotates.append(alpha);
            rotAngles.phiRotates.append(azimut);
        }
        initFrame = frame;
        file.close();
    }
    else
    {
        throw std::runtime_error("Не удалось открыть файл " + filename.toStdString());
    }
}


void MLSTask::readLines(const QString& filename, bool skipFirstRow, Catalog& catalog)
{
    lineFrame.clear();
    QFile file (filename);
    if (file.open(QIODevice::ReadOnly))
    {
        QTextStream in(&file);
        QString line;
        constexpr const quint16 rowLength = 8;
        if (skipFirstRow)
        {
            in.readLineInto(&line);
        }
        QString previousAlpha;
        QString previousAzimut;

        in.readLineInto(&line);
        QStringList splitted = line.split("\t");
        previousAzimut = splitted[2];
        LineData lineData;
        auto fillLineData = [&catalog, &lineData, this](double alpha, double delta, QStringList splitted)
        {
            lineData.x.append((splitted[5].toDouble() - frameX) * pixelSize);
            lineData.y.append((splitted[6].toDouble() - frameX) * pixelSize);
            ++lineData.count;
            for (int i = 0; i < catalog.alpha.size(); i++)
            {
                if (qFuzzyCompare(catalog.alpha[i], alpha)
                        && qFuzzyCompare(catalog.delta[i], delta))
                {
                    lineData.vecSt.append(i);
                }
            }

        };
        fillLineData(splitted[3].toDouble(), splitted[4].toDouble(), splitted);

        while (in.readLineInto(&line))
        {
            QStringList list = line.split("\t");
            if (list.size() != rowLength)
            {
                throw std::logic_error("Неверная строка в каталоге");
            }
            if (previousAzimut != list[2])
            {
                lineFrame.append(lineData);
                lineData.clear();
                previousAzimut = list[2];
            }
            fillLineData(list[3].toDouble(), list[4].toDouble(), list);
        }
        lineFrame.append(lineData);
        file.close();
    }
    else
    {
        throw std::runtime_error("Не удалось открыть файл " + filename.toStdString());
    }
}


void MLSTask::readLinesWithCatalog(const QString& filename, bool skipFirstRow, Catalog& catalog)
{
    QFile file (filename);
    if (file.open(QIODevice::ReadOnly))
    {
        QTextStream in(&file);
        QString line;
        //constexpr const quint16 rowLength = 8;
        if (skipFirstRow)
        {
            in.readLineInto(&line);
        }
        qint32 pos = in.pos();
        QString previousAlpha;
        QString previousAzimut;

        in.readLineInto(&line);
        QStringList splitted = line.split("\t");
        previousAzimut = splitted[2];
        in.seek(pos);
        bool finish = false;
        double l, m, n;
        while (!finish)
        {
            in.readLineInto(&line);
            splitted = line.split("\t", QString::SkipEmptyParts);
            if (!splitted.isEmpty() && previousAzimut == splitted[2])
            {
                catalog.alpha.append(splitted[3].toDouble());
                catalog.delta.append(splitted[4].toDouble());
                calculateDirectionVector(catalog.alpha.last(), catalog.delta.last(), l, m, n);
                catalog.L.append(l);
                catalog.M.append(m);
                catalog.N.append(n);
            }
            else
            {
                finish = true;
            }
        }
    }
    readLines(filename, skipFirstRow, catalog);
}

void MLSTask::calculateXY(double Mstand[3][3], double collimator[3], double focus, double& X, double& Y)
{
    double LMN[3];
    LMN[0] = Mstand[0][0] * collimator[0] + Mstand[0][1] * collimator[1] + Mstand[0][2] * collimator[2];
    LMN[1] = Mstand[1][0] * collimator[0] + Mstand[1][1] * collimator[1] + Mstand[1][2] * collimator[2];
    LMN[2] = Mstand[2][0] * collimator[0] + Mstand[2][1] * collimator[1] + Mstand[2][2] * collimator[2];
    if (!fishEye)
    {
        X = -focus * LMN[0] / LMN[2];
        Y = -focus * LMN[1] / LMN[2];
        //qDebug() << X << Y << "nofisheye";
    }
    else
    {
        X = (focus * atan2m(sqrt(pow(LMN[0], 2) + pow(LMN[1], 2)), LMN[2]))
                / sqrt(pow(LMN[1] / LMN[0], 2) + 1);
        Y = (focus * atan2m(sqrt(pow(LMN[0], 2) + pow(LMN[1], 2)), LMN[2]))
                / sqrt(pow(LMN[0] / LMN[1], 2) + 1);
        //qDebug() << X << Y << "fisheye";
    }

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
        case REVERSE_X:
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
        case REVERSE_Y:
            m[0][0] = 1;
            m[1][1] = -1;
            m[2][2] = 1;
            break;
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
void MLSTask::calculateFull(const QBitArray& derivativeFlags, Results& results, ResultErrors& errors)
{
    calculatePrivate(derivativeFlags, results, errors, rotAngles, frame);

}

void MLSTask::firstApprox(const Catalog& catalog, const LineData& lineData, OrientAngles& angles)
{
    qint32 i ,j, imin, imax, jmin, vmin, vmax;
    double min, max, dist, min_dist, CA, SA, CX;
    double cos_al[2], sin_al[2], cos_dlt[2], sin_dlt[2], cos_st, sin_st;

    min = 1000.; max = -1000.;  min_dist = 1000.;

    for (int j = 0; j < lineData.count; j++)
    {
        //star with min and max delta
        i = lineData.vecSt[j];
        if (catalog.N[i] > max)  {max = catalog.N[i]; imax = j;}
        if (catalog.N[i] < min)  {min = catalog.N[i]; imin = j;}
        //nearest for the center
        dist = sqrtm(lineData.x[j] * lineData.x[j] + lineData.y[j] * lineData.y[j]);
        if (dist < min_dist) { min_dist = dist; jmin = j; }
    }
    if (imax == imin)
    {
        imax = lineData.count - 1;
    }

    //в качестве начального приближения альфа и дельта - параметры ближайшей к центру звезды
    j = lineData.vecSt[jmin];
    angles.dl = asinm(catalog.N[j]);

    angles.al = atan2m(catalog.M[j],catalog.L[j]);
    if (angles.al < 0) angles.al += 2 * M_PI;

    //первое приближение азимута - по двум крайним звездам
    vmin = lineData.vecSt[imin];
    vmax = lineData.vecSt[imax];
    sin_dlt[0] = catalog.N[vmin];
    sin_dlt[1] = catalog.N[vmax];

    cos_dlt[0] = 1 - sin_dlt[0] * sin_dlt[0]; cos_dlt[0] = sqrtm(cos_dlt[0]);
    cos_dlt[1] = 1 - sin_dlt[1] * sin_dlt[1]; cos_dlt[1] = sqrtm(cos_dlt[1]);

    sin_al[0] = catalog.M[vmin] / cos_dlt[0];
    cos_al[0] = catalog.L[vmin] / cos_dlt[0];

    sin_al[1] = catalog.M[vmax] / cos_dlt[1];
    cos_al[1] = catalog.L[vmax] / cos_dlt[1];

    cos_st = catalog.L[vmin] * catalog.L[vmax]
            + catalog.M[vmin] * catalog.M[vmax]
            + catalog.N[vmin] * catalog.N[vmax];
    sin_st = 1 - cos_st * cos_st; sin_st = sqrtm(sin_st);

    CA = (sin_dlt[0] - sin_dlt[1] * cos_st) / (sin_st * cos_dlt[1]);
    SA = cos_dlt[0] * (sin_al[0] * cos_al[1] - cos_al[0] * sin_al[1]) / sin_st;

    angles.Az = atan2m(SA, CA);
    if (angles.Az < 0) angles.Az += 2 * M_PI;
    SA = lineData.x[imin] - lineData.x[imax]; CA = lineData.y[imin] - lineData.y[imax];

    CX = atan2m(SA,CA); if (CX < 0) CX += 2 * M_PI;
    angles.Az = CX - angles.Az; if (angles.Az < 0) angles.Az += 2 * M_PI;
    angles.Az += M_PI;
}

void MLSTask::calculateOldModelPrivate(const Catalog& catalog, const LineData& lineData, OldErrors errors, double& focus, qint32 numP, bool dist, bool save)
{
    focus = 32;
    int clk = 0;
    double mxy=1000;
    double mxy_pr, mx, my;
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
    double CC1, CC2, CC3, CC4, CC5;
    double cos_al, sin_al, cos_dlt, sin_dlt, cos_Az, sin_Az;
    OrientAngles orientAngles;
    firstApprox(catalog, lineData, orientAngles);
    double Mst[3][3];
    do
    {
        mxy_pr = mxy;
        mxy = 0.;  mx = 0.; my = 0.;
        calcTransitionMatrix(orientAngles.al * radToDegrees,
                             orientAngles.dl * radToDegrees,
                             orientAngles.Az * radToDegrees,
                             Mst);
        cos_al=cos(orientAngles.al);   sin_al=sin(orientAngles.al);
        cos_dlt=cos(orientAngles.dl);  sin_dlt=sin(orientAngles.dl);
        cos_Az=cos(orientAngles.Az);   sin_Az=sin(orientAngles.Az);

        for (int j = 0; j < lineData.count; j++)
        {
            int i = lineData.vecSt[j];
            int k = j * 2; // номер уравнения для х
            int l = k + 1; // номер уравнения для у

            CC1 = Mst[2][0] * catalog.L[i] + Mst[2][1] * catalog.M[i] + Mst[2][2] * catalog.N[i];
            CC2 = Mst[2][1] * catalog.L[i] - Mst[2][0] * catalog.M[i];
            CC3 = (catalog.L[i] * cos_al + catalog.M[i] * sin_al) * sin_dlt - catalog.N[i] * cos_dlt;
            CC4 = -focus * ((Mst[0][0] * catalog.L[i] + Mst[0][1] * catalog.M[i] + Mst[0][2] * catalog.N[i])/CC1);
            CC5 = -focus * ((Mst[1][0] * catalog.L[i] + Mst[1][1] * catalog.M[i] + Mst[1][2] * catalog.N[i])/CC1);

            DFR[k] = (double)(lineData.x[j] - CC4);
            DFR[l] = (double)(lineData.y[j] - CC5);

            if (numP == 1)
            {
                DRV[k][0] = (double)(CC4 / focus);
                DRV[l][0] = (double)(CC5 / focus);
            }
            else
            {
                DRV[k][0] = (double)(CC5);
                DRV[l][0] = (double)(-CC4);
                DRV[k][1] = (double)(focus*sin_Az+CC4*(CC3 / CC1));
                DRV[l][1] = (double)(focus*cos_Az+CC5*(CC3 / CC1));
                DRV[k][2] = (double)((focus*(Mst[0][1]*catalog.L[i]-Mst[0][0]*catalog.M[i])+CC4*CC2)/CC1);
                DRV[l][2] = (double)((focus*(Mst[1][1]*catalog.L[i]-Mst[1][0]*catalog.M[i])+CC5*CC2)/CC1);
            }
            if (numP > 3)
            {
                DRV[k][3]=(double)(CC4 / focus);
                DRV[l][3]=(double)(CC5 / focus);
            }

            mxy += DFR[k] * DFR[k] + DFR[l] * DFR[l];
            mx  += DFR[k] * DFR[k];
            my  += DFR[l] * DFR[l];
        }

        mxy = mxy / (double)(2 * lineData.count - numP);    mxy = sqrtm(mxy);
        mx  = mx / (double)(lineData.count - numP / 2.);    mx  = sqrtm(mx);
        my  = my / (double)(lineData.count - numP / 2.);    my  = sqrtm(my);

        for (int i = 0; i<numP; i++)
        {
            for (int j = 0; j<numP; j++) // произведение обычной и транс. матрицы параметров
            {
                DRVM[i][j] = 0.;      //(DRV)^T*DRV
                for (int k = 0; k < 2 * lineData.count; k++)
                {
                    DRVM[i][j] += DRV[k][i] * DRV[k][j];
                }
            }
        }

        //        for (int i = 0; i < numP; i++)
        //        {
        //            QString line;
        //            for (int j = 0; j < numP; j++)
        //            {
        //                line.append(QString::number(DRVM[i][j]) + " ");
        //            }
        //            qDebug() << line;
        //        }

        if (numP == 1)
            DRVM_1[0][0] = 1/DRVM[0][0];
        else
            invertMatrix(DRVM, DRVM_1, numP); // вычисление обратной матрицы


        for (int i = 0; i < numP; i++) { // умножение транспонированной на вектор расхождений координат
            DRVH[i] = 0.;             //DRVH=(DRV)^T*DIF
            for (int k = 0; k < 2 * lineData.count; k++) DRVH[i] += DRV[k][i] * DFR[k];
        }

        for (int i = 0; i < numP; i++)
        { // вычисление дельт углов
            ANGD[i] = 0.;
            for (int k = 0; k < numP; k++)
            {
                ANGD[i] += DRVM_1[i][k] * DRVH[k];
            }
        }

        if (numP==1) focus += ANGD[0];
        else
        {
            orientAngles.Az+=ANGD[0];
            orientAngles.dl+=ANGD[1];
            orientAngles.al+=ANGD[2];
        }

        if (numP>3) focus+=ANGD[3];
        clk++;
    }
    while ((fabs(mxy - mxy_pr) > maxError) && (clk < maxIterations));

    orientAngles.dl  = asin(Mst[2][2]);
    orientAngles.al  = atan2m(Mst[2][1], Mst[2][0]);    if (orientAngles.al < 0)  orientAngles.al += 2 * M_PI;
    orientAngles.Az  = atan2m(Mst[0][2], Mst[1][2]);    if (orientAngles.Az < 0) orientAngles.Az += 2 * M_PI;

    errors.dAl  = sqrtm(DRVM_1[2][2]) * mxy * radToSec * cos(orientAngles.dl);
    errors.dDlt = sqrtm(DRVM_1[1][1]) * mxy * radToSec;
    errors.dAz  = sqrtm(DRVM_1[0][0]) * mxy * radToSec;

    if (numP > 3)
    {
        errors.dFoc = sqrtm(DRVM_1[3][3]) * mxy * 1000.;
    }
    else errors.dFoc = 0.;
    qDebug() <<"Orient:" << orientAngles.al * radToDegrees << orientAngles.dl * radToDegrees << orientAngles.Az * radToDegrees;
    qDebug() <<"Error:" << errors.dFoc <<  errors.dAl  << errors.dDlt <<  errors.dAz <<  mxy * 1000 << mx * 1000  << my * 1000 << "\n";

    if (dist)
    {
        for (int i = 0; i < lineData.count; i++)
        {
            int deltaStep = i * 2;
            distData.x.append(lineData.x[i]);
            distData.y.append(lineData.y[i]);
            distData.dx.append(DFR[deltaStep]);
            distData.dy.append(DFR[deltaStep + 1]);
        }
    }
    if (save)
    {
        QFile experimentalData("experimental_data.txt");
        if (experimentalData.open(QIODevice::Append))
        {
            QTextStream out(&experimentalData);
            if (experimentalData.size() == 0)
            {
                out << QString("foc\talpha\tdelta\tazimut\terAlpha\terDelta\terAzimut\tmx\tmy\tmxy\tmodel\n");
            }
            out << QString("%1\t%2\t%3\t%4\t%5\t%6\t%7\t%8\t%9\t%10\t%11\n")
                   .arg(focus)
                   .arg(orientAngles.al * radToDegrees).arg(orientAngles.dl * radToDegrees).arg(orientAngles.Az * radToDegrees)
                   .arg(errors.dAl).arg(errors.dDlt).arg(errors.dAz)
                   .arg(mx * 1000).arg(my * 1000).arg(mxy * 1000)
                   .arg(modelName);
        }

    }
    for (int count = 0; count < maxDer; count++)
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

void MLSTask::calculateOldModel(const Catalog& catalog, OldErrors errors, double focus, qint32 numP)
{
    double initFocus = focus;
    QVector <double> focuses;
    for (const auto& i : lineFrame)
    {
        calculateOldModelPrivate(catalog, i, errors, focus, numP, false, false);
        focuses.append(focus);
        focus = initFocus;
    }
    focus = calculateMean(focuses.begin(), focuses.end(), 0.0);
    for (const auto& i : lineFrame)
    {
        calculateOldModelPrivate(catalog, i, errors, focus, numP - 1, true, false);
    }
    findDistorsio(3, true);
    for (int i = 0; i < lineFrame.size(); i++)
    {
        for (int j = 0; j < lineFrame[i].count; j++)
        {
            lineFrame[i].x[j] = calculateDistorsio(lineFrame[i].x[j], lineFrame[i].x[j], lineFrame[i].y[j], xDistV);
            lineFrame[i].y[j] = calculateDistorsio(lineFrame[i].y[j], lineFrame[i].x[j], lineFrame[i].y[j], yDistV);
        }
    }
    for (const auto& i : lineFrame)
    {
        calculateOldModelPrivate(catalog, i, errors, focus, numP - 1, false, true);
    }
    saveDistorsio();
}


void MLSTask::readStandCatalog(const QString& filename, Catalog& catalog)
{
    QFile file (filename);
    if(file.open(QIODevice::ReadOnly))
    {
        QTextStream in(&file);
        QString line;
        constexpr const quint16 rowLength = 5;
        double l, m, n;
        while(in.readLineInto(&line))
        {
            QStringList list = line.split("\t");
            if (list.size() != rowLength)
                throw std::logic_error("Неверная строка в каталоге");
            catalog.alpha.append(list[3].toDouble());
            catalog.delta.append(list[4].toDouble());
            calculateDirectionVector(catalog.alpha.last(), catalog.delta.last(), l, m, n);
            catalog.L.append(l);
            catalog.M.append(m);
            catalog.N.append(n);
        }
        catalog.countStars = catalog.alpha.size();
        file.close();
    }
}

void MLSTask::calculateTriangle(bool distorsio, Results& results, ResultErrors& errors)
{
    constexpr const double angleMainFirst = 0.353553;
    constexpr const double angleMainSecond = 0.353553;
    constexpr const double angleFirstSecond = 0.5;
    constexpr const double angleBetwTriangles = 2.0;
    const QVector <double> fixedAngles {angleMainFirst, angleMainSecond, angleFirstSecond, angleBetwTriangles};

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
    double mxy = 1000;
    double mxy_pr;
    qint32 clk = 0;
    qint32 numP = 0;
    double focus = results.foc;
    QList <double> distX;
    QList <double> distY;
    QList <double> distXDrv;
    QList <double> distYDrv;
    for (int i = 0; i < 21; i++)
    {
        distX.append(0);
        distY.append(0);
    }
    distXDrv = distX;
    distYDrv = distY;
    qint32 counter = 0;
    qint32 nCft = 10;
    //double distDelta = 0.0001;
    QVector <double > distDelta {1.0e-12, 1.0e-10, 1.0e-10, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6};
    //попробовать сделать 20 коэф (10 - х, 10 - y) и попробовать уровняться со всеми возможными комбинациями коэффициентов
    do
    {
        mxy_pr = mxy;
        mxy = 0.;
        counter = 0;
        for (const auto& line : triangleFrame)
        {

            for (int i = 0; i < line.size(); i++)
            {
                for (int j = 0; j < line[i].points.size() - 1; j++)
                {
                    for (int k = j + 1; k < line[i].points.size(); k++)
                    {
                        if (!distorsio)
                        {
                            double angle = calculateAngle(line[i].points[j], line[i].points[k], focus, distX, distY);
                            DFR[counter] = fixedAngles[j + k - 1] - angle;
                            double nAngle = calculateAngle(line[i].points[j], line[i].points[k], focus + deltaFocus, distX, distY);
                            DRV[counter][numP] = (nAngle - angle) / deltaFocus;
                            mxy += DFR[counter] * DFR[counter];
                            ++counter;
                        }
                        else
                        {
                            numP = 0;
                            for (int l = 0; l < nCft; l++)
                            {
                                distXDrv = distX;
                                distXDrv[l] += distDelta[l];
                                double angle = calculateAngle(line[i].points[j], line[i].points[k], focus, distX, distY);
                                DFR[counter] = fixedAngles[j + k - 1] - angle;
                                double nAngle = calculateAngle(line[i].points[j], line[i].points[k], focus, distXDrv, distY);
                                DRV[counter][numP++] = (nAngle - angle) / distDelta[l];
                                mxy += DFR[counter] * DFR[counter];
                                distYDrv = distY;
                                distYDrv[l] += distDelta[l];
                                DFR[counter] = angle;
                                nAngle = calculateAngle(line[i].points[j], line[i].points[k], focus, distX, distYDrv);
                                DRV[counter][numP++] = (nAngle - angle) / distDelta[l];
                                mxy += DFR[counter] * DFR[counter];
                            }
                            ++counter;
                        }

                    }
                }
                for (int j = i + 1; j < line.size(); j++)
                {
                    if (i == line.size() - 1)
                    {
                        break;
                    }
                    if (!distorsio)
                    {
                        double angle = abs(calculateAngle(line[i].points[0], line[j].points[0], focus, distX, distY));
                        DFR[counter] = fixedAngles.last() * (j - i) - angle;
                        double nAngle = calculateAngle(line[i].points[0], line[j].points[0], focus + deltaFocus, distX, distY);
                        DRV[counter][numP] = (abs(nAngle) - angle) / deltaFocus;
                        mxy += DFR[counter] * DFR[counter];
                        ++counter;
                    }
                    else
                    {
                        numP = 0;
                        for (int l = 0; l < nCft; l++)
                        {
                            distXDrv = distX;
                            distXDrv[l] += distDelta[l];
                            double angle = abs(calculateAngle(line[i].points[0], line[j].points[0], focus, distX, distY));
                            DFR[counter] = fixedAngles.last() * (j - i) - angle;;
                            double nAngle = abs(calculateAngle(line[i].points[0], line[j].points[0], focus, distXDrv, distY));
                            DRV[counter][numP++] = (nAngle - angle) / distDelta[l];
                            mxy += DFR[counter] * DFR[counter];

                            distYDrv = distY;
                            distYDrv[l] += distDelta[l];
                            DFR[counter] = angle;
                            nAngle = calculateAngle(line[i].points[0], line[j].points[0], focus, distX, distYDrv);
                            DRV[counter][numP++] = (nAngle - angle) / distDelta[l];
                            mxy += DFR[counter] * DFR[counter];
                        }
                        ++counter;
                    }
                }
            }
        }
        if (!distorsio)
        {
            numP = 1;
        }
        else
        {
            //                        for (int i = 0; i < counter; i++)
            //                        {
            //                            QString str;
            //                            for (int j = 0; j < numP; j++)
            //                            {
            //                                str.append(QString::number(DRV[i][j]) + " ");
            //                            }
            //                            qDebug() << str;
            //                        }
        }

        mxy = mxy / (double)(counter - numP);  mxy = sqrtm(mxy);
        for (int i = 0; i < numP; i++)
        {
            for (int j = 0; j < numP; j++) // произведение обычной и транс. матрицы параметров
            {
                DRVM[i][j] = 0.;      //(DRV)^T*DRV
                for (int k = 0; k < counter; k++)
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
            invertMatrix(DRVM, DRVM_1, numP); // вычисление обратной матрицы
        }
        for (int i = 0; i < numP; i++)
        { // умножение транспонированной на вектор расхождений координат
            DRVH[i] = 0.;             //DRVH=(DRV)^T*DIF
            for (int k = 0; k < counter; k++)
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

        if (!distorsio)
        {
            focus += ANGD[0];
        }
        else
        {
            for (int i = 0; i < nCft; i++)
            {
                distX[i] += ANGD[2 * i];
                distY[i] += ANGD[2 * i + 1];
            }
        }
        clk++;
    }
    while ((fabs(mxy - mxy_pr) > maxError)
           && (clk < maxIterations));

    if (!distorsio)
    {
        results.foc = focus;
        errors.dfoc = sqrtm(DRVM_1[0][0]) * mxy;
        qDebug() << results.foc << errors.dfoc << mxy;
    }
    else
    {
        qDebug() << distX << distY;
    }



    if (mxy > thershold * pixelSize)
    {
        throw std::runtime_error("Грубые измерения.");
    }

    for (int count = 0; count < maxDer; count++)
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


    double collimator[3] {-(cos(ang.lambdaOY) * cos(ang.lambdaOX)),
                -(cos(ang.lambdaOY) * sin(ang.lambdaOX)),
                -sin(ang.lambdaOY)};
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
            invertMatrix(DRVM, DRVM_1, numP); // вычисление обратной матрицы
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

    for (int count = 0; count < maxDer; count++)
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


void MLSTask::findDistorsio(int nPow, bool clear)
{
    findDistCft(nPow, distData.x, distData.y, distData.dx, distData.dy);
    if (clear)
    {
        distData.x.clear();
        distData.y.clear();
        distData.dx.clear();
        distData.dy.clear();
    }
}


void MLSTask::saveDistorsio()
{
    QFile file(QString("dist%1.txt").arg(modelName));
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
        out << QString("x_init y_init dx dy\n");
        for (int i = 0; i < distData.x.size(); i++)
        {
            out << QString("%1\t%2\t%3\t%4\n")
                   .arg(initFrame[i].x())
                   .arg(initFrame[i].y())
                   .arg(distData.dx[i])
                   .arg(distData.dy[i]);
        }
    }
}

void MLSTask::includeDistorsio()
{

    for (int i = 0; i < frame.size(); i++)
    {
        auto tmpFrame = frame[i];
        frame[i].setX(calculateDistorsio(tmpFrame.x(), tmpFrame.x(), tmpFrame.y(), xDistV));
        frame[i].setY(calculateDistorsio(tmpFrame.y(), tmpFrame.x(), tmpFrame.y(), yDistV));
    }
}

QString MLSTask::printTestTable(const QString& filename, bool dist, double focus)
{
    double frameXF = frameX * 2;
    double frameYF = frameY * 2;
    QVector <QPointF> points {
                QPointF(0, 0),
                QPointF(frameXF, 0),
                QPointF(frameXF, frameYF),
                QPointF(0, frameYF),
                QPointF(0, frameYF / 2),
                QPointF(frameXF, frameYF / 2),
                QPointF(frameXF / 2, 0),
                QPointF(frameXF / 2, frameYF)
                /*,
                        QPointF(frameXF * (1./2), frameYF * (1./3)),
                        QPointF(frameXF * (2./3), frameYF * (1./2)),
                        QPointF(frameXF * (1./2), frameYF * (2./3)),
                        QPointF(frameXF * (1./3), frameYF * (1./2))*/};

    for(int i = 0; i < points.size(); i++)
    {
        points[i].setX((points[i].x() - frameX) * pixelSize);
        points[i].setY((points[i].y() - frameY) * pixelSize);
    }
    if (dist)
    {
        for (int i = 0; i < points.size(); i++)
        {
            QPointF tmpPoint = points[i];
            points[i].setX(calculateDistorsio(tmpPoint.x(), tmpPoint.x(), tmpPoint.y(), xDistV));
            points[i].setY(calculateDistorsio(tmpPoint.y(), tmpPoint.x(), tmpPoint.y(), yDistV));
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
        }
    }
    QString retLine;
    retLine.append("Диагональ сверху-вниз: " + QString::number(m[0][2]) + "\n");
    retLine.append("Диагональ снизу-вверх: " + QString::number(m[1][3]) + "\n");
    retLine.append("Левая сторона: " + QString::number(m[0][1]) + "\n");
    retLine.append("Правая сторона: " + QString::number(m[2][3]) + "\n");
    retLine.append("Верхняя сторона: " + QString::number(m[0][3]) + "\n");
    retLine.append("Нижняя сторона: " + QString::number(m[2][1]) + "\n");
    retLine.append("Нижняя сторона: " + QString::number(m[2][1]) + "\n");
    retLine.append("Центр сверху-вниз: " + QString::number(m[4][5]) + "\n");
    retLine.append("Центр слева-направо: " + QString::number(m[6][7]) + "\n");
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
    return retLine;
}

ShiftData MLSTask::getShiftData() const
{
    ShiftData data;
    for (int i = 0; i < initFrame.size(); i++)
    {
        data.x.append(initFrame[i].x());
        data.y.append(initFrame[i].y());
        data.dx.append(distData.dx[i]);
        data.dy.append(distData.dy[i]);
        data.azimut.append(rotAngles.phiRotates[i]);
    }
    return data;
}

void MLSTask::findDistCft(int Npow, QVector <double>& x, QVector <double>& y, QVector <double>& dx, QVector <double>& dy)
{
    double CX[maxParams], CY[maxParams],B[maxParams][maxParams],
            B_1[maxParams][maxParams], AX[maxParams], AY[maxParams];
    double Mcft_i[maxParams];
    double x_2, y_2, x_3, y_3, x_4, y_4, x_5,
            y_5, x_6, y_6, x_7, y_7, x_8, y_8, x_9, y_9;
    int Ncft;

    switch (Npow)
    {
        case 2: Ncft = 6;  break;
        case 3: Ncft = 10; break;
        case 4: Ncft = 15; break;
        case 5: Ncft = 21; break;
        case 6: Ncft = 28; break;
        case 7: Ncft = 36; break;
        case 8: Ncft = 45; break;
        case 9: Ncft = 55; break;
        default: Npow = 3; Ncft = 10; break;
    }


    for (int i = 0; i < maxParams; i++)
    {
        CX[i] = 0.; CY[i] = 0.;
        for (int j = 0; j < Ncft; j++)
            B[i][j] = 0.;
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
            x_3 = x_2 * x[i];
            y_3 = y_2 * y[i];

            Mcft_i[6] = x_3;
            Mcft_i[7] = x_2 * y[i];
            Mcft_i[8] = x[i] * y_2;
            Mcft_i[9] = y_3;
            if (Npow >= 4)
            {
                x_4 = x_3 * x[i];
                y_4 = y_3 * y[i];

                Mcft_i[10] = x_4;
                Mcft_i[11] = x_3 * y[i];
                Mcft_i[12] = x_2 * y_2;
                Mcft_i[13] = x[i] * y_3;
                Mcft_i[14] = y_4;

                if (Npow >= 5)
                {
                    x_5 = x_4 * x[i];
                    y_5 = y_4 * y[i];

                    Mcft_i[15] = x_5;
                    Mcft_i[16] = x_4 * y[i];
                    Mcft_i[17] = x_3 * y_2;
                    Mcft_i[18] = x_2 * y_3;
                    Mcft_i[19] = x[i] * y_4;
                    Mcft_i[20] = y_5;
                    if (Npow >= 6)
                    {
                        x_6 = x_5 * x[i];
                        y_6 = y_5 * y[i];

                        Mcft_i[21] = x_6;
                        Mcft_i[22] = x_5 * y[i];
                        Mcft_i[23] = x_4 * y_2;
                        Mcft_i[24] = x_3 * y_3;
                        Mcft_i[25] = x_2 * y_4;
                        Mcft_i[26] = x[i] * y_5;
                        Mcft_i[27] = y_6;
                        if (Npow >= 7)
                        {
                            x_7 = x_6 * x[i];
                            y_7 = y_6 * y[i];

                            Mcft_i[28] = x_7;
                            Mcft_i[29] = x_6 * y[i];
                            Mcft_i[30] = x_5 * y_2;
                            Mcft_i[31] = x_4 * y_3;
                            Mcft_i[32] = x_3 * y_4;
                            Mcft_i[33] = x_2 * y_5;
                            Mcft_i[34] = x[i] * y_6;
                            Mcft_i[35] = y_7;
                            if (Npow >= 8)
                            {
                                x_8=x_7 * x[i];
                                y_8=y_7 * y[i];

                                Mcft_i[36] = x_8;
                                Mcft_i[37] = x_7 * y[i];
                                Mcft_i[38] = x_6 * y_2;
                                Mcft_i[39] = x_5 * y_3;
                                Mcft_i[40] = x_4 * y_4;
                                Mcft_i[41] = x_3 * y_5;
                                Mcft_i[42] = x_2 * y_6;
                                Mcft_i[43] = x[i] * y_7;
                                Mcft_i[44] = y_8;
                                if (Npow >= 9)
                                {
                                    x_9 = x_8 * x[i];
                                    y_9 = y_8 * y[i];

                                    Mcft_i[45] = x_9;
                                    Mcft_i[46] = x_8 * y[i];
                                    Mcft_i[47] = x_7 * y_2;
                                    Mcft_i[48] = x_6 * y_3;
                                    Mcft_i[49] = x_5 * y_4;
                                    Mcft_i[50] = x_4 * y_5;
                                    Mcft_i[51] = x_3 * y_6;
                                    Mcft_i[52] = x_2 * y_7;
                                    Mcft_i[53] = x[i] * y_8;
                                    Mcft_i[54] = y_9;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (int k = 0; k < Ncft; k++)
        {
            CX[k] += Mcft_i[k] * dx[i];
            CY[k] += Mcft_i[k] * dy[i];
            for (int j = 0; j < Ncft; j++)
                B[k][j] += Mcft_i[k] * Mcft_i[j];
        }

    }

    gaussObr(Ncft,B,B_1);    //(A^T*A)^(-1)

    for (int i = Ncft; i < maxParams; i++)
    {
        AX[i] = 0.; AY[i] = 0;
    }

    for (int i = 0; i < Ncft; i++)
    {
        AX[i] = 0.; AY[i] = 0.;
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

    distData.dxDiff.clear();
    distData.dyDiff.clear();
    for (int i = 0; i < x.size(); i++)
    {

        x_2 = x[i]*x[i];
        y_2= y[i]*y[i];
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
        distData.dxDiff.append(dx_new);
        distData.dyDiff.append(dy_new);
    }
}


int MLSTask::gaussObr(int cntStar, double mass[55][55], double mObr[55][55])
{
    int i,j,k;
    double a,b;
    double sum;

    for (i = 0;i < cntStar; i++)
    {
        for(j = 0;j < cntStar; j++)
            mObr[i][j] = 0;
        mObr[i][i] = 1;
    }

    for (i = 0;i < cntStar; i++)
    {
        a = mass[i][i];
        if (fabs(a) < 1e-38) return 0;
        for (j= i+1; j < cntStar; j++)
        {
            b=mass[j][i];
            for (k = 0; k < cntStar; k++)
            {
                mass[j][k] = mass[j][k] - mass[i][k] * b / a;
                mObr[j][k] = mObr[j][k] - mObr[i][k] * b / a;
            }
        }
    }

    for (i = 0; i < cntStar; i++)
    {
        for (j = cntStar - 1; j >= 0; j--)
        {
            sum=0;
            for (k = cntStar - 1; k > j; k--)
                sum += mass[j][k] * mObr[k][i];
            if (qFuzzyCompare(mass[j][j], 0))
            {
                return 0;
            }
            mObr[j][i] = (mObr[j][i] - sum) / mass[j][j];
        }
    }
    return 1;
}

double MLSTask::calculateAngle(const QPointF& fPoint, const QPointF& sPoint, double focus, QList <double>& distorsioCoefX, QList <double>& distorsioCoefY)
{
    QPointF fPointDist (fPoint.x() + calculateDistorsioDelta(fPoint.x(), fPoint.y(), distorsioCoefX),
                        fPoint.y() + calculateDistorsioDelta(fPoint.x(), fPoint.y(), distorsioCoefY));
    double fLength = sqrtm (pow(fPointDist.x(), 2) + pow(fPointDist.y(), 2)  + pow(focus, 2));
    fPointDist.setX(fPointDist.x() / fLength);
    fPointDist.setY(fPointDist.y() / fLength);
    double normFocusf = -focus / fLength;
    QPointF sPointDist (sPoint.x() + calculateDistorsioDelta(sPoint.x(), sPoint.y(), distorsioCoefX),
                        sPoint.y() + calculateDistorsioDelta(sPoint.x(), sPoint.y(), distorsioCoefY));
    double sLength = sqrtm (pow(sPointDist.x(), 2) + pow(sPointDist.y(), 2)  + pow(focus, 2));
    sPointDist.setX(sPointDist.x() / sLength);
    sPointDist.setY(sPointDist.y() / sLength);
    double normFocuss = -focus / sLength;
    double cos = calculateScalarProduct(fPointDist.x(), sPointDist.x(), fPointDist.y(), sPointDist.y(), normFocusf, normFocuss);
    return acosm(cos) * radToDegrees;
}
