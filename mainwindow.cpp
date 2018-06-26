#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    settings = new QSettings(QDir::currentPath() + "/" + "config.ini", QSettings::IniFormat, this);
    loadSettings();
    ui->powComboBox->addItems(QStringList {"3", "5"});
    ui->formatTypeComboBox->addItems(QStringList{"Формат №1", "Формат №2(Никитин)"});
    //    double quat1[4]{-0.419788873899061,
    //                 -0.238363459945942,
    //                 0.649585552399954,
    //                 0.587365840358774};
    //    double angles[3];

    //     quatToEkvA(quat1,angles);
    //     qDebug() << angles[0] << angles[1] << angles[2];
    //    double quat2[4]{-0.419564124451634,
    //                 -0.236331216961869,
    //                 0.649326466321183,
    //                 0.588632582859084
    //};
    //        double angles2[3];
    //         quatToEkvA(quat2,angles2);
    //         qDebug() << angles2[0] << angles2[1] << angles2[2];
    //            qDebug() << angles[0] - angles2[0] << angles[1] -angles2[1] << angles[2] -angles2[2];
    //    double quat1[4]{0.499440916946064,
    //                    0.650944714273445,
    //                    -0.098704654996217,
    //                    -0.563104687639594
    //};
    //    double angles[3];

    //     quatToEkvA(quat1,angles);
    //     qDebug() << angles[0] << angles[1] << angles[2];
    //    double quat2[4]{0.499414246077221,
    //                    0.6481685866302,
    //                    -0.0978453471311587,
    //                    -0.566470825829316

    //};
    //        double angles2[3];
    //         quatToEkvA(quat2,angles2);
    //         qDebug() << angles2[0] << angles2[1] << angles2[2];
    //            qDebug() << angles[0] - angles2[0] << angles[1] -angles2[1] << angles[2] -angles2[2];
    //    double quat1[4]{0.212915879772215,
    //                    0.657417678099037,
    //                    -0.0375695997619592,
    //                    0.721842941833932

    //};
    //    double angles[3];

    //     quatToEkvA(quat1,angles);
    //     qDebug() << angles[0] << angles[1] << angles[2];
    //    double quat2[4]{0.209797806226135,
    //                    0.657261057540879,
    //                    -0.0391216354117838,
    //                    0.722815404801227
    //};
    //        double angles2[3];
    //         quatToEkvA(quat2,angles2);
    //         qDebug() << angles2[0] << angles2[1] << angles2[2];
    //            qDebug() << angles[0] - angles2[0] << angles[1] -angles2[1] << angles[2] -angles2[2];



    //    double quat1[4]{    			0.2110919506,
    //                            -0.0562356822,
    //                                -0.9637190108,
    //                                -0.1533735246};
    //double angles[3];

    //quatToEkvA(quat1,angles);
    //qDebug() << angles[0] << angles[1] << angles[2];
    //double quat2[4]{0.2110987903,
    //                -0.0562792868,
    //                -0.9637168650,
    //                -0.1533615143
    //               };
    //double angles2[3];
    //quatToEkvA(quat2,angles2);
    //qDebug() << angles2[0] << angles2[1] << angles2[2];
    //qDebug() << angles[0] - angles2[0] << angles[1] -angles2[1] << angles[2] -angles2[2];
}

MainWindow::~MainWindow()
{
    saveSettings();
    delete ui;
}

void MainWindow::saveSettings()
{
    settings->setValue("treshold", ui->measureThesholdSpinBox->value());
    settings->setValue("xMatrix", ui->xMatrixSpinBox->value());
    settings->setValue("yMatrix", ui->yMatrixSpinBox->value());
    settings->setValue("pixSize", ui->pixSizeSpinBox->value());
    settings->setValue("focCheck", ui->focusCheckBox->isChecked());
    settings->setValue("focValue", ui->focusSpinBox->value());
    settings->setValue("lambdaYCheck", ui->lambdaYCheckBox->isChecked());
    settings->setValue("lambdaYValue", ui->lambdaOYSpinBox->value());
    settings->setValue("lambdaXCheck", ui->lambdaXCheckBox->isChecked());
    settings->setValue("lambdaXValue", ui->lambdaOXSpinBox->value());
    settings->setValue("gammaXCheck", ui->gammaXCheckBox->isChecked());
    settings->setValue("gammaXValue", ui->gammaOXSpinBox->value());
    settings->setValue("gammaYCheck", ui->gammaYCheckBox->isChecked());
    settings->setValue("gammaYValue", ui->gammaOYSpinBox->value());
    settings->setValue("gammaZCheck", ui->gammaZCheckBox->isChecked());
    settings->setValue("gammaZValue", ui->gammaOZSpinBox->value());
    settings->setValue("alpha2YCheck", ui->alpha2YCheckBox->isChecked());
    settings->setValue("alpha2YValue", ui->alpha2OYSpinBox->value());
    settings->setValue("alpha2XCheck", ui->alpha2XCheckBox->isChecked());
    settings->setValue("alpha2XValue", ui->alpha2OXSpinBox->value());
    settings->setValue("path",ui->pathLineEdit->text());

    settings->sync();

}

void MainWindow::loadSettings()
{
    ui->measureThesholdSpinBox->setValue(settings->value("treshold").toInt());
    ui->xMatrixSpinBox->setValue(settings->value("xMatrix").toInt());
    ui->yMatrixSpinBox->setValue(settings->value("yMatrix").toInt());
    ui->pixSizeSpinBox->setValue(settings->value("pixSize").toDouble());
    ui->focusCheckBox->setChecked(settings->value("focCheck").toBool());
    ui->focusSpinBox->setValue(settings->value("focValue").toDouble());
    ui->lambdaYCheckBox->setChecked(settings->value("lambdaYCheck").toBool());
    ui->lambdaOYSpinBox->setValue(settings->value("lambdaYValue").toDouble());
    ui->lambdaXCheckBox->setChecked(settings->value("lambdaXCheck").toBool());
    ui->lambdaOXSpinBox->setValue(settings->value("lambdaXValue").toDouble());
    ui->gammaXCheckBox->setChecked(settings->value("gammaXCheck").toBool());
    ui->gammaOXSpinBox->setValue(settings->value("gammaXValue").toDouble());
    ui->gammaYCheckBox->setChecked(settings->value("gammaYCheck").toBool());
    ui->gammaOYSpinBox->setValue(settings->value("gammaYValue").toDouble());
    ui->gammaZCheckBox->setChecked(settings->value("gammaZCheck").toBool());
    ui->gammaOZSpinBox->setValue(settings->value("gammaZValue").toDouble());
    ui->alpha2YCheckBox->setChecked(settings->value("alpha2YCheck").toBool());
    ui->alpha2OYSpinBox->setValue(settings->value("alpha2YValue").toDouble());
    ui->alpha2XCheckBox->setChecked(settings->value("alpha2XCheck").toBool());
    ui->alpha2OXSpinBox->setValue(settings->value("alpha2XValue").toDouble());
    ui->pathLineEdit->setText(settings->value("path").toString());
}

void MainWindow::on_pushButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Open modeling file"),
                                                    QDir::currentPath(), tr("*.txt"));
    QDir::setCurrent(fileName);
    ui->pathLineEdit->setText(fileName);
}
Results MainWindow::setFirstApprox()
{
    Results results;
    results.foc = ui->focusSpinBox->value();
    results.lambdaOY = ui->lambdaOYSpinBox->value();
    results.lambdaOX = ui->lambdaOXSpinBox->value();
    results.alphaTwoOY = ui->alpha2OYSpinBox->value();
    results.alphaTwoOX = ui->alpha2OXSpinBox->value();
    results.gammaOZ = ui->gammaOZSpinBox->value();
    results.gammaOY = ui->gammaOYSpinBox->value();
    results.gammaOX = ui->gammaOXSpinBox->value();
    return results;
}

QBitArray MainWindow::setFlags()
{
    QBitArray derivativeFlags(DERIVATIVES::COUNT);
    derivativeFlags.setBit(DERIVATIVES::LAMBDA_OY, ui->lambdaYCheckBox->isChecked());
    derivativeFlags.setBit(DERIVATIVES::LAMBDA_OX, ui->lambdaXCheckBox->isChecked());
    derivativeFlags.setBit(DERIVATIVES::ALPHA_TWO_OY, ui->alpha2YCheckBox->isChecked());
    derivativeFlags.setBit(DERIVATIVES::ALPHA_TWO_OX, ui->alpha2XCheckBox->isChecked());
    derivativeFlags.setBit(DERIVATIVES::GAMMA_OZ, ui->gammaZCheckBox->isChecked());
    derivativeFlags.setBit(DERIVATIVES::GAMMA_OY, ui->gammaYCheckBox->isChecked());
    derivativeFlags.setBit(DERIVATIVES::GAMMA_OX, ui->gammaXCheckBox->isChecked());
    derivativeFlags.setBit(DERIVATIVES::PHI, false);
    derivativeFlags.setBit(DERIVATIVES::FOCUS, ui->focusCheckBox->isChecked());
    return derivativeFlags;
}


void MainWindow::printResults(Results& res)
{
    ui->textEdit->append("\n=================================================");
    ui->textEdit->append("\nРезультаты: \n");
    ui->textEdit->append(QString("%1%2%3%4%5%6%7%8")
                         .arg("Фокус", -10, ' ').arg("ЛямбдаОУ", 10, ' ').arg("ЛямбдаОХ", 10, ' ')
                         .arg("Альфа2ОУ", 10, ' ').arg("Альфа2ОХ", 10, ' ')
                         .arg("ГаммаOZ", 10, ' ').arg("ГаммаОУ", 10, ' ').arg("ГаммаОХ", 10, ' '));
    ui->textEdit->append(QString("%1%2%3%4%5%6%7%8")
                         .arg(res.foc, -10, 'f', 4).arg(res.lambdaOY, 10, 'f', 4)
                         .arg(res.lambdaOX, 10, 'f', 4).arg(res.alphaTwoOY, 10, 'f', 4)
                         .arg(res.alphaTwoOX, 10, 'f', 4).arg(res.gammaOZ, 10, 'f', 4)
                         .arg(res.gammaOY, 10, 'f', 4)
                         .arg(res.gammaOX, 10, 'f', 4));

}
void MainWindow::printErrors(ResultErrors& err)
{
    ui->textEdit->append("\nОшибки (углы в секундах, остальное в микронах): \n ");
    ui->textEdit->append(QString("%1%2%3%4%5%6%7%8%9%10%11")
                         .arg("Фокус", 10, ' ').arg("ЛямбдаОУ", 10, ' ').arg("ЛямбдаОХ", 10, ' ')
                         .arg("Альфа2ОУ", 10, ' ').arg("Альфа2ОХ", 10, ' ')
                         .arg("ГаммаOZ", 10, ' ').arg("ГаммаОУ", 10, ' ').arg("ГаммаОХ", 10, ' ')
                         .arg("mx", 10, ' ').arg("my", 10, ' ').arg("mxy", 10, ' '));
    ui->textEdit->append(QString("%1%2%3%4%5%6%7%8%9%10%11")//%1|%2|%3|%4|%5|%6|%7|%8|%9|%10|%11
                         .arg(err.dfoc, 10, 'f', 2, ' ')
                         .arg(err.dlambdaOY, 10, 'f', 2, ' ').arg(err.dlambdaOX, 12, 'f', 2, ' ')
                         .arg(err.dalphaTwoOY, 12, 'f', 2, ' ').arg(err.dalphaTwoOX, 14, 'f', 2, ' ')
                         .arg(err.dgammaOZ, 14, 'f', 2, ' ').arg(err.dgammaOY, 14, 'f', 2, ' ').arg(err.dgammaOX, 13, 'f', 2, ' ')
                         .arg(err.mx, 10, 'f', 2, ' ').arg(err.my, 10, 'f', 2, ' ').arg(err.mxy, 10, 'f', 2, ' '));
    ui->textEdit->append("\n=================================================");

}

void MainWindow::printAngles (const QString& before, const QString& after)
{
    ui->textEdit->append("\nУглы до дисторсии:\n");
    ui->textEdit->append(before);
    if (!after.isEmpty())
    {
        ui->textEdit->append("\nУглы после дисторсии:\n");
        ui->textEdit->append(after);
    }

}
void MainWindow::saveErrors(ResultErrors& err)
{
    QFile file("errors.txt");
    if (file.open(QIODevice::WriteOnly))
    {
        QTextStream out(&file);
        out << QString("%1%2%3%4%5%6%7%8%9%10%11")
               .arg("Фокус", -10, ' ').arg("ЛямбдаОУ", -10, ' ').arg("ЛямбдаОХ",-10, ' ')
               .arg("Альфа2ОУ", -10, ' ').arg("Альфа2ОХ", -10, ' ')
               .arg("ГаммаОZ", -10, ' ').arg("ГаммаОУ", -10, ' ').arg("ГаммаОХ", -10, ' ')
               .arg("mx", -10, ' ').arg("my", -10, ' ').arg("mxy", -10, ' ') << "\n";

        out << QString("%1%2%3%4%5%6%7%8%9%10%11")//%1|%2|%3|%4|%5|%6|%7|%8|%9|%10|%11
               .arg(err.dfoc, 10, 'f', 2, ' ')
               .arg(err.dlambdaOY, 10, 'f', 2, ' ').arg(err.dlambdaOX, 12, 'f', 2, ' ')
               .arg(err.dalphaTwoOY, 12, 'f', 2, ' ').arg(err.dalphaTwoOX, 14, 'f', 2, ' ')
               .arg(err.dgammaOZ, 14, 'f', 2, ' ').arg(err.dgammaOY, 14, 'f', 2, ' ').arg(err.dgammaOX, 13, 'f', 2, ' ')
               .arg(err.mx, 10, 'f', 2, ' ').arg(err.my, 10, 'f', 2, ' ').arg(err.mxy, 10, 'f', 2, ' ');

    }
}
void MainWindow::saveResults(Results& res)
{
    QFile file("results.txt");
    if (file.open(QIODevice::WriteOnly))
    {
        QTextStream out(&file);
        out << QString("%1%2%3%4%5%6%7%8")
               .arg("Фокус", 10, ' ').arg("ЛямбдаОУ", 10, ' ').arg("ЛямбдаОХ", 10, ' ')
               .arg("Альфа2ОУ", 10, ' ').arg("Альфа2ОХ", 10, ' ')
               .arg("ГаммаОУ", 10, ' ').arg("ГаммаОХ", 10, ' ') << "\n";
        out << QString("%1%2%3%4%5%6%7%8")
               .arg(res.foc, 8, ' ').arg(res.lambdaOY, 8, ' ')
               .arg(res.lambdaOX, 8, ' ').arg(res.alphaTwoOY, 8, ' ')
               .arg(res.alphaTwoOX, 8, ' ').arg(res.gammaOZ, 8, 'f')
               .arg(res.gammaOY, 8, ' ')
               .arg(res.gammaOX, 8, ' ');;

    }
}

bool MainWindow::checkChoose()
{
    return  ui->lambdaYCheckBox->isChecked()
            || ui->lambdaXCheckBox->isChecked()
            || ui->alpha2YCheckBox->isChecked()
            || ui->alpha2XCheckBox->isChecked()
            || ui->gammaYCheckBox->isChecked()
            || ui->gammaXCheckBox->isChecked()
            || ui->focusCheckBox->isChecked();
}
void MainWindow::on_pushButton_2_clicked()
{
    try
    {
        if (!checkChoose())
        {
            ui->textEdit->append("Не выбраны параметры для расчёта\n");
            return;
        }

        quint32 frameX = ui->xMatrixSpinBox->value() / 2;
        quint32 frameY = ui->yMatrixSpinBox->value() / 2;
        double pixSize = ui->pixSizeSpinBox->value();
        task.clearAll();
        task.setFrameSize(frameX, frameY);
        task.setPixelSize(pixSize);
        task.setMeasureTheshold(ui->measureThesholdSpinBox->value());
        if (ui->modelRadioButton->isChecked())
        {
            task.readModelData(ui->pathLineEdit->text(), true);
        }
        else
        {
            task.readRealData(ui->pathLineEdit->text(), true);
        }

        if (ui->upRadioButton->isChecked())
        {
            task.setAxisDirection(UP);
        }
        else if (ui->downRadioButton->isChecked())
        {
            task.setAxisDirection(DOWN);
        }
        else if (ui->leftRadioButton->isChecked())
        {
            task.setAxisDirection(LEFT);
        }
        else if (ui->rightRadioButton->isChecked())
        {
            task.setAxisDirection(RIGHT);
        }

        auto flags = setFlags();
        auto results = setFirstApprox();
        ResultErrors errors;

        if (ui->fitFocusCheckBox->isChecked()
                && flags.at(DERIVATIVES::FOCUS))
        {
            task.fitFocusByLines(flags, results, errors);
            flags.setBit(DERIVATIVES::FOCUS, false);
        }

        task.calculate(flags, results, errors);
        task.saveShifts("before_dist");
        auto angBeforeDist = task.printTestTable("ang_distance_before_dist.txt", false, results.foc);
        QVector <QString> angAfterDist;
        if (ui->distRadioButton->isChecked())
        {
            quint32 pow = ui->powComboBox->currentText().toInt();
            task.findDistorsio(pow);
            task.saveDistorsio();
            task.includeDistorsio();
            flags.setBit(DERIVATIVES::FOCUS, false);
            task.calculate(flags, results, errors);
            angAfterDist = task.printTestTable("ang_distance_after_dist.txt", true, results.foc);
            ui->textEdit->append("\nУглы после дисторсии:\n");
            ui->textEdit->append(angAfterDist[0]);
            task.saveShifts("after_dist");
        }
        printResults(results);
        printErrors(errors);

        if (angAfterDist.size() > 0)
        {
            printAngles(angBeforeDist[0], angAfterDist[0]);
        }
        else
        {
            printAngles(angBeforeDist[0], QString());
        }
        saveResults(results);
        saveErrors(errors);
    }
    catch (std::exception& e)
    {
        ui->textEdit->append(QString(e.what()) + "\n");
    }

}

void MainWindow::on_pushButton_3_clicked()
{
    ui->textEdit->clear();
}


void MainWindow::joinMeasureFiles()
{

}

void MainWindow::on_chooseRawFilesPushButton_clicked()
{
    QStringList files = QFileDialog::getOpenFileNames(
                this,
                "Select one or more files to open",
                "/home",
                "Text files (*.txt)");
    QStringList coordFiles;
    QStringList angleFiles;

    copy_if(files.constBegin(), files.constEnd(), back_inserter(angleFiles),
            [](auto& i){return i.contains("protocols");});
    copy_if(files.constBegin(), files.constEnd(), back_inserter(coordFiles),
            [](auto& i){return i.contains("REZ");});
    if (coordFiles.size() == 0 || angleFiles.size() == 0)
    {
        ui->textEdit->append("Нечего объединять. Один из типов файлов отсутствует\n");
        return;
    }

    //    auto parse = [this](auto& v, const auto& pos) {
    //        for (const auto& i : v)
    //        {
    //            QFile file(i);
    //            if (file.open(QIODevice::ReadOnly))
    //            {
    //                QTextStream in(&file);
    //                QString line;
    //                while (in.readLineInto(&line))
    //                {
    //                    QStringList list = line.split("\t", QString::SplitBehavior::SkipEmptyParts);
    //                    v.append(list);
    //                }
    //            }
    //            else
    //            {
    //                ui->textEdit->append("Не удалось открыть файл" + i + "\n");
    //                return;
    //            }
    //            file.close();
    //        }
    //        for (int i = 0; i < v.size(); i++)
    //        {

    //            quint32 pos = v[i][pos].indexOf(QRegExp("\\\\(?:(?!\\\\))\\d\\d\\."));
    //            v[i][pos] = v[i][pos].mid(pos);
    //        }
    //    };
    QVector <QStringList> coordFilesData;
    // parse(coordFilesData, 13);
    for (const auto& i : coordFiles)
    {
        QFile file(i);
        if (file.open(QIODevice::ReadOnly))
        {
            QTextStream in(&file);
            QString line;
            while (in.readLineInto(&line))
            {
                QStringList list = line.split("\t", QString::SplitBehavior::SkipEmptyParts);
                coordFilesData.append(list);
            }
        }
        else
        {
            ui->textEdit->append("Не удалось открыть файл" + i + "\n");
            return;
        }
        file.close();
    }
    for (int i = 0; i < coordFilesData.size(); i++)
    {

        quint32 pos =
                coordFilesData[i][13].indexOf(QRegExp("\\\\(?:(?!\\\\))\\d\\d\\."));
        coordFilesData[i][13] =
                coordFilesData[i][13].mid(pos);
    }

    QVector <QStringList> anglesFilesData;
    //parse(anglesFilesData, 1);
    for (const auto&i : angleFiles)
    {
        QFile file(i);
        if (file.open(QIODevice::ReadOnly))
        {
            QTextStream in(&file);
            QString line;
            while (in.readLineInto(&line))
            {
                QStringList list = line.split("\t", QString::SplitBehavior::SkipEmptyParts);
                anglesFilesData.append(list);
            }
        }
        else
        {
            ui->textEdit->append("Не удалось открыть файл" + i + "\n");
            return;
        }
        file.close();
    }
    for (int i = 0; i < anglesFilesData.size(); i++)
    {

        quint32 pos =
                anglesFilesData[i][1].indexOf(QRegExp("\\\\(?:(?!\\\\))\\d\\d\\."));
        anglesFilesData[i][1] =
                anglesFilesData[i][1].mid(pos);
    }
    sort(coordFilesData.begin(), coordFilesData.end(), [](auto& a, auto& b){return a[13] < b[13];});
    sort(anglesFilesData.begin(), anglesFilesData.end(), [](auto& a, auto& b){return a[1] < b[1];});

    QVector <QStringList> resFilesData;
    quint32 pos = 0;
    for (int i = 0; i < coordFilesData.size(); i++)
    {
        for (int j = pos; j < anglesFilesData.size(); j++)
        {
            if (coordFilesData[i][13] == anglesFilesData[j][1])
            {
                anglesFilesData[j][2].replace(",",".");
                anglesFilesData[j][3].replace(",",".");
                resFilesData.append(coordFilesData[i]);
                resFilesData.last().append(anglesFilesData[j][2]);
                resFilesData.last().append(anglesFilesData[j][3]);
                pos = j;
                break;
            }
        }
    }

    pos = files[0].indexOf(QRegExp("(/)(?!.+/)"), 0);
    QDir dir(files[0].mid(0, pos));
    QFile file(dir.absolutePath() + "/" + QDateTime::currentDateTime().toString("dd.MM.yyyy hh-mm-ss") + "result_file_union.txt");
    if (file.open(QIODevice::WriteOnly))
    {
        QTextStream out(&file);
        if (ui->formatTypeComboBox->currentIndex() == 0)
        {
            out << "X\t" << "Y\t" << "Alpha\t" << "Azimut\n";
            for (const auto& i : resFilesData)
            {
                out << QString("%1\t%2\t%3\t%4\n")
                       .arg(i[1]).arg(i[2])
                        .arg(i[14]).arg(i[15]);
            }
        }
        else
        {
            out << "Data	Time	File	X	Y	I	N	maxI	Aldeg	Gm deg\n";
            for (const auto& i : resFilesData)
            {
                out << QString("%1\t%2\t%3\t%4\t%5\t%6\t%7\t%8\t%9\t%10\n")
                       .arg("123").arg("123").arg("123").arg(i[1]).arg(i[2])
                        .arg("123").arg("123").arg("123")
                        .arg(i[14]).arg(i[15]);
            }
        }
    }
}

