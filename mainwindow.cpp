#include "mainwindow.h"
#include "ui_mainwindow.h"



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent), ui(new Ui::MainWindow),
    deleteSc(QKeySequence("Delete"), this)
{
    ui->setupUi(this);
    connect(&deleteSc, &QShortcut::activated, this, &MainWindow::on_removeToolButton_clicked);
    settings = new QSettings(QDir::currentPath() + "/" + "config.ini", QSettings::IniFormat, this);
    loadSettings();
    ui->powComboBox->addItems(QStringList {"2", "3", "5"});
    ui->formatTypeComboBox->addItems(QStringList{"Формат №1", "Формат №2(Никитин)"});
    plotter.reset (new CXYPlotter(ui->plot));
    plotter->setGraphName("Рассогласования");
    plotter->setAxisXName("X, пикс");
    plotter->setAxisYName("Y, пикс");
    plotter->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectItems | QCP::iSelectPlottables| QCP::iSelectLegend);
    plotter->setSelectable(QCP::SelectionType::stNone);
    plotter->setUseLegend(false);
    plotter->setAutoReplot(false);
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
    settings->setValue("reverse", ui->reverseCheckBox->isChecked());
    settings->setValue("real_data", ui->realDataRadioButton->isChecked());
    settings->setValue("model_data", ui->modelRadioButton->isChecked());
    settings->setValue("scale", ui->scaleSpinBox->value());
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
    ui->reverseCheckBox->setChecked(settings->value("reverse").toBool());
    ui->realDataRadioButton->setChecked(settings->value("real_data").toBool());
    ui->modelRadioButton->setChecked(settings->value("model_data").toBool());
    ui->scaleSpinBox->setValue(settings->value("scale").toInt());
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
    ui->textEdit->append("<font color=green>\nРезультаты: \n</font>");
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
    ui->textEdit->append("\n<font color=green>Ошибки (углы в секундах, остальное в микронах):</font> \n ");
    ui->textEdit->append(QString("%1%2%3%4%5%6%7%8%9%10%11")
                         .arg("Фокус").arg("ЛямбдаОУ", 10, ' ').arg("ЛямбдаОХ", 10, ' ')
                         .arg("Альфа2ОУ", 10, ' ').arg("Альфа2ОХ", 10, ' ')
                         .arg("ГаммаOZ", 10, ' ').arg("ГаммаОУ", 10, ' ').arg("ГаммаОХ", 10, ' ')
                         .arg("mx", 10, ' ').arg("my", 10, ' ').arg("mxy", 10, ' '));
    ui->textEdit->append(QString("%1%2%3%4%5%6%7%8%9%10%11")//%1|%2|%3|%4|%5|%6|%7|%8|%9|%10|%11
                         .arg(err.dfoc, 0, 'f', 2)
                         .arg(err.dlambdaOY, 10, 'f', 2, ' ').arg(err.dlambdaOX, 12, 'f', 2, ' ')
                         .arg(err.dalphaTwoOY, 12, 'f', 2, ' ').arg(err.dalphaTwoOX, 14, 'f', 2, ' ')
                         .arg(err.dgammaOZ, 14, 'f', 2, ' ').arg(err.dgammaOY, 14, 'f', 2, ' ').arg(err.dgammaOX, 13, 'f', 2, ' ')
                         .arg(err.mx, 10, 'f', 2, ' ').arg(err.my, 10, 'f', 2, ' ').arg(err.mxy, 10, 'f', 2, ' '));
    ui->textEdit->append("\n=================================================");

}

void MainWindow::printAngles (const QString& before, const QString& after)
{
    ui->textEdit->append("\n<font color=green>Углы до дисторсии:</font>\n");
    ui->textEdit->append(before);
    if (!after.isEmpty())
    {
        ui->textEdit->append("\n<font color=green>Углы после дисторсии:</font>\n");
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
               .arg(res.gammaOX, 8, ' ');

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
            [](auto& i) {return i.contains("protocols");});
    copy_if(files.constBegin(), files.constEnd(), back_inserter(coordFiles),
            [](auto& i) {return i.contains("REZ");});
    if (coordFiles.size() == 0 || angleFiles.size() == 0)
    {
        ui->textEdit->append("Нечего объединять. Один из типов файлов отсутствует\n");
        return;
    }
    QVector <QStringList> coordFilesData;
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
        qint32 pos =
                coordFilesData[i][13].indexOf(QRegExp("(\\\\)(?!.*\\\\)"));
        coordFilesData[i][13] =
                coordFilesData[i][13].mid(pos);
    }

    QVector <QStringList> anglesFilesData;
    for (const auto&i : angleFiles)
    {
        QFile file(i);
        if (file.open(QIODevice::ReadOnly))
        {
            QTextStream in(&file);
            QString line;
            qint32 count = 0;
            while (in.readLineInto(&line))
            {
                QStringList list = line.split("\t", QString::SplitBehavior::SkipEmptyParts);
                anglesFilesData.append(list);
                anglesFilesData.last().append("\t" + QString::number(count++));
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

        qint32 pos =
                anglesFilesData[i][1].indexOf(QRegExp("(\\\\)(?!.*\\\\)"));
        anglesFilesData[i][1] =
                anglesFilesData[i][1].mid(pos);
    }
    sort(coordFilesData.begin(), coordFilesData.end(), [](auto& a, auto& b){return a[13] < b[13];});
    sort(anglesFilesData.begin(), anglesFilesData.end(), [](auto& a, auto& b){return a[1] < b[1];});
    QVector <QStringList> resFilesData;
    qint32 pos = files[0].indexOf(QRegExp("(/)(?!.+/)"), 0);
    QDir dir(files[0].mid(0, pos));
    bool flag = false;
    if (ui->threeFileCheckBox->isChecked())
    {
        for (int i = 0; i < coordFilesData.size(); i++)
        {
            for (int j = 0; j < anglesFilesData.size(); j++)
            {
                if (coordFilesData[i][13] == anglesFilesData[j][1])
                {
                    flag = true;
                    anglesFilesData[j][2].replace(",", ".");
                    anglesFilesData[j][3].replace(",", ".");
                    resFilesData.append(coordFilesData[i]);
                    resFilesData.last().append(anglesFilesData[j][2]);
                    resFilesData.last().append(anglesFilesData[j][3]);
                    resFilesData.last().append(anglesFilesData[j].last());
                    qDebug() << resFilesData.last();
                }
            }
            if (flag == false)
            {
                qDebug() << coordFilesData[i][13] ;
            }
            flag = false;
        }


        QFile file(dir.absolutePath() + "/" + QDateTime::currentDateTime().toString("dd.MM.yyyy hh-mm-ss") + "result_file_union_three.txt");
        if (file.open(QIODevice::WriteOnly))
        {
            QTextStream out(&file);
            if (ui->formatTypeComboBox->currentIndex() == 0)
            {
                out << "X\t" << "Y\t" << "Alpha\t" << "Azimut\t" << "Cluster №\n";
                for (const auto& i : resFilesData)
                {
                    out << QString("%1\t%2\t%3\t%4\t%5\n")
                           .arg(i[1])
                            .arg(i[2])
                            .arg(i[14])
                            .arg(i[15])
                            .arg(i.last());
                }
            }
        }
    }
    else
     {
        for (int i = 0; i < coordFilesData.size(); i++)
        {
            for (int j = 0; j < anglesFilesData.size(); j++)
            {
                if (coordFilesData[i][13] == anglesFilesData[j][1])
                {
                    flag = true;
                    anglesFilesData[j][2].replace(",", ".");
                    anglesFilesData[j][3].replace(",", ".");
                    resFilesData.append(coordFilesData[i]);
                    resFilesData.last().append(anglesFilesData[j][2]);
                    resFilesData.last().append(anglesFilesData[j][3]);
                    break;
                }
            }
            if (flag == false)
            {
                qDebug() << coordFilesData[i][13] ;
            }
            flag = false;
        }


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

}



void MainWindow::on_saveToolButton_clicked()
{
    if (editStarted)
    {
        QFile file ("modifiedPointList.txt");
        if (file.open(QIODevice::WriteOnly))
        {
            QTextStream out (&file);
            out << "X	Y	Alpha	Azimut\n";
            for (int i = 0; i < editingList.size(); i++)
            {
                out << editingList[i] << "\n";
            }
        }
        file.close();
    }
    editStarted = false;
}

void MainWindow::on_removeToolButton_clicked()
{
    if (editStarted)
    {
        auto list = ui->plot->selectedItems();
        if (!list.isEmpty())
        {
            Q_ASSERT(selectedIndex != -1);
            QCPItemTracer* tracer = static_cast <QCPItemTracer*> (list.first());
            ui->plot->removeGraph(tracer->graph());
            ui->plot->removeItem(tracer);
            editingList.removeAt(selectedIndex);
            selectedIndex = -1;
            ui->plot->replot();
        }
    }
    else
    {
        QMessageBox::information(nullptr, tr("Ошибка"), tr("Сначала начните редактирование исходного файла."));
    }
}

void MainWindow::on_startToolButton_clicked()
{
    QFile file (ui->pathLineEdit->text());
    if (file.open(QIODevice::ReadOnly))
    {
        QTextStream in (&file);
        QString line;
        in.readLineInto(&line);
        while (in.readLineInto(&line))
        {
            editingList.append(line);
        }
        file.close();
        editStarted = true;
    }
    else
    {
        QMessageBox::information(nullptr, tr("Ошибка"), tr("Не удалось открыть исходный файл."));
    }
}

void MainWindow::on_chooseCatData_clicked()
{

    Catalog catalog;
    qint32 frameX = ui->xMatrixSpinBox->value() / 2;
    qint32 frameY = ui->yMatrixSpinBox->value() / 2;
    double pixSize = ui->pixSizeSpinBox->value();
    task.clearAll();
    task.setFrameSize(frameX, frameY);
    task.setPixelSize(pixSize);
    task.setMeasureTheshold(ui->measureThesholdSpinBox->value());

    QStringList files = QFileDialog::getOpenFileNames(
                this,
                "Select one or more files to open",
                "",
                "TXT (*.txt)");
    for (const auto& i : files)
    {
        OldErrors err;
        qint32 stPos = i.lastIndexOf("/") + 1;
        QString modelName = i.mid(stPos, i.lastIndexOf(".txt") - stPos);
        task.setModelName(modelName);
        Catalog catalog;
        task.readLinesWithCatalog(i, true, catalog);
        task.calculateOldModel(catalog, err, ui->focusSpinBox->value(), 4);
    }
}


void MainWindow::on_startCalcPushButton_clicked()
{
    try
    {
        if (!checkChoose())
        {
            ui->textEdit->append("Не выбраны параметры для расчёта\n");
            return;
        }
        qint32 frameX = ui->xMatrixSpinBox->value() / 2;
        qint32 frameY = ui->yMatrixSpinBox->value() / 2;
        double pixSize = ui->pixSizeSpinBox->value();
        task.clearAll();
        task.setFrameSize(frameX, frameY);
        task.setPixelSize(pixSize);
        task.setMeasureTheshold(ui->measureThesholdSpinBox->value());
        if (ui->modelRadioButton->isChecked())
        {
            if (ui->allFrameRadioButton->isChecked())
            {
                task.readFullModelData(ui->pathLineEdit->text(), true);
            }
            else
            {
                task.readTriangleModelData(ui->pathLineEdit->text(), true);
            }
        }
        else
        {
            task.readRealData(ui->pathLineEdit->text(), true, ui->reverseCheckBox->isChecked());
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
        else if (ui->mirrorYRadioButton->isChecked())
        {
            task.setAxisDirection(REVERSE_Y);
        }
        else if (ui->mirrorXRadioButton->isChecked())
        {
            task.setAxisDirection(REVERSE_X);
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

        if (ui->allFrameRadioButton->isChecked())
        {
            task.calculateFull(flags, results, errors);
        }
        else
        {
            task.calculateTriangle(false, results, errors);
            task.calculateTriangle(true, results, errors);
        }
        task.saveShifts("before_dist");
        auto angBeforeDist = task.printTestTable("ang_distance_before_dist.txt", false, results.foc);
        QString angAfterDist;
        if (ui->distRadioButton->isChecked())
        {
            qint32 pow = ui->powComboBox->currentText().toInt();
            task.findDistorsio(pow);
            task.saveDistorsio();
            task.includeDistorsio();
            qint32 setFlags = 0;
            for (int i = 0; i < flags.count(); i++)
            {
                if (flags.at(i))
                {
                    setFlags++;
                }
            }
            if (setFlags == 1 && flags.at(DERIVATIVES::FOCUS))
            {
                flags.setBit(DERIVATIVES::FOCUS, true);
            }
            else
            {
                flags.setBit(DERIVATIVES::FOCUS, false);
            }
            task.calculateFull(flags, results, errors);
            angAfterDist = task.printTestTable("ang_distance_after_dist.txt", true, results.foc);
            ui->textEdit->append("\nУглы после дисторсии:\n");
            task.saveShifts("after_dist");
        }
        printResults(results);
        printErrors(errors);

        if (angAfterDist.size() > 0)
        {
            printAngles(angBeforeDist, angAfterDist);
        }
        else
        {
            printAngles(angBeforeDist, QString());
        }
        saveResults(results);
        saveErrors(errors);

        ShiftData shiftData = task.getShiftData();
        ui->plot->clearPlottables();
        ui->plot->clearItems();
        plotter->setRangeAxisX(-frameX, frameX);
        plotter->setRangeAxisY(-frameY, frameY);
        plotter->setStyle(QCPGraph::LineStyle::lsLine, QCPScatterStyle::ssStar);

        for (int i = 0; i < shiftData.dx.size(); i++)
        {
            shiftData.x[i] = shiftData.x[i] / pixSize; //+ frameX;
            shiftData.dx[i] = shiftData.dx[i] / pixSize; //+ frameX;
            shiftData.dx[i] *= ui->scaleSpinBox->value();
            shiftData.y[i] = shiftData.y[i] / pixSize;// + frameY;
            shiftData.dy[i] = shiftData.dy[i] / pixSize;// + frameY;
            shiftData.dy[i] *= ui->scaleSpinBox->value();
            QVector <double> x {shiftData.x[i], shiftData.x[i] + shiftData.dx[i]};
            QVector <double> y {shiftData.y[i], shiftData.y[i] + shiftData.dy[i]};
            plotter->addDefaultGraph(y, x, false);
            qint32 colorIndex = static_cast <qint32> (std::abs(shiftData.azimut[i]));
            QCPItemTracer* tracer = plotter->setTracer(ui->plot->graphCount() - 1, x[0], colorIndex, true);
            connect(tracer, &QCPItemTracer::selectionChanged, [this, tracer](bool selected) {
                if (selected)
                {
                    for (int i = 0; i < ui->plot->graphCount(); i++)
                    {
                        if (ui->plot->graph(i) == tracer->graph())
                        {
                            selectedIndex = i;
                            ui->choosedPointInfoLabel->setText(QString("X: %1 Y: %2 Индекс: %3")
                                                               .arg(tracer->graphKey() + (ui->xMatrixSpinBox->value() / 2))
                                                               .arg(tracer->graph()->dataMainValue(0) + (ui->yMatrixSpinBox->value() / 2))
                                                               .arg(selectedIndex));
                        }
                    }
                }
            });
        }
        ui->plot->replot();


    }
    catch (std::exception& e)
    {
        ui->textEdit->append(QString(e.what()) + "\n");
    }

}

void MainWindow::on_clearTextEditPushButton_clicked()
{
    ui->textEdit->clear();
}

void MainWindow::on_chooseModelFilePushButton_clicked()
{
    static QString fn = QString();
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Open modeling file"),
                                                    fn, tr("*.txt"));
    fn = fileName;
    QString dirPath = fileName;
    qint32 pos = dirPath.indexOf(QRegExp("(/)(?!.*/)"));
    QDir dir(dirPath.remove(pos, fileName.size() - pos));
    dir.setCurrent(dir.absolutePath());
    ui->pathLineEdit->setText(fileName);
}










//    //QFile file1 ("//Camera20/data/БОКЗ-М60-1000/БОКЗ-М60-1000 №12/Небо три прибора - 2018.10.31/30.10.2018 22.50.46 - 12 13 14/quatinfo_f1.txt");
//    QFile file1 ("//Camera20/data/БОКЗ-М60-1000/БОКЗ-М60-1000 №12/Небо три прибора - 2018.11.16/16.11.2018 21.12.18 - 12 13 14 повороты15 угл.мин-с + вибро/quatinfo_1.txt");
//    QVector <QPair <int, QVector<double>>> f;
//    if (file1.open(QIODevice::ReadOnly))
//    {
//        QTextStream in(&file1);
//        QString line;
//        while (in.readLineInto(&line))
//        {
//            QStringList list = line.split("\t");

//            f.append(qMakePair <int, QVector<double>>(list[0].toInt(), QVector<double> {list[1].toDouble(), list[2].toDouble(), list[3].toDouble(), list[4].toDouble(),
//                                                                                        list[5].toDouble(), list[6].toDouble()/*, list[7].toDouble(), list[8].toDouble(), list[9].toDouble()*/}));
//        }
//    }
//    //QFile file2 ("//Camera20/data/БОКЗ-М60-1000/БОКЗ-М60-1000 №12/Небо три прибора - 2018.10.31/30.10.2018 22.50.46 - 12 13 14/quatinfo_s2.txt");
//    QFile file2 ("//Camera20/data/БОКЗ-М60-1000/БОКЗ-М60-1000 №12/Небо три прибора - 2018.11.16/16.11.2018 21.12.18 - 12 13 14 повороты15 угл.мин-с + вибро/quatinfo_2.txt");
//    QVector <QPair <int, QVector<double>>> s;
//    if (file2.open(QIODevice::ReadOnly))
//    {
//        QTextStream in(&file2);
//        QString line;
//        while (in.readLineInto(&line))
//        {
//            QStringList list = line.split("\t");
//            s.append(qMakePair<int, QVector<double>>(list[0].toInt(), QVector<double>{list[1].toDouble(), list[2].toDouble(), list[3].toDouble(), list[4].toDouble(),
//                                                                                      list[5].toDouble(), list[6].toDouble()/*, list[7].toDouble(), list[8].toDouble(), list[9].toDouble()*/}));
//        }
//    }
//    //QFile file3 ("//Camera20/data/БОКЗ-М60-1000/БОКЗ-М60-1000 №12/Небо три прибора - 2018.10.31/30.10.2018 22.50.46 - 12 13 14/quatinfo_t3.txt");
//    QFile file3 ("//Camera20/data/БОКЗ-М60-1000/БОКЗ-М60-1000 №12/Небо три прибора - 2018.11.16/16.11.2018 21.12.18 - 12 13 14 повороты15 угл.мин-с + вибро/quatinfo_3.txt");
//    QVector <QPair <int, QVector<double>>> t;
//    if (file3.open(QIODevice::ReadOnly))
//    {
//        QTextStream in(&file3);
//        QString line;
//        while (in.readLineInto(&line))
//        {
//            QStringList list = line.split("\t");
//            t.append(qMakePair<int, QVector<double>>(list[0].toInt(), QVector<double> {list[1].toDouble(), list[2].toDouble(), list[3].toDouble(), list[4].toDouble(),
//                                                                                       list[5].toDouble(), list[6].toDouble()/*, list[7].toDouble(), list[8].toDouble(), list[9].toDouble()*/}));
//        }
//    }
//    std::sort(f.begin(), f.end(), [](auto& a, auto& b) {return a.first < b.first;});
//    std::sort(s.begin(), s.end(), [](auto& a, auto& b) {return a.first < b.first;});
//    std::sort(t.begin(), t.end(), [](auto& a, auto& b) {return a.first < b.first;});

//    QFile outf("out90.txt");
//    if (outf.open(QIODevice::WriteOnly))
//    {
//        QTextStream out(&outf);
//        out << "TimePr\t3-1 OZ\t3-2 OZ\t2-1 OZ\t3-1 XY\t3-2 XY\t2-1 XX\tdet1\tloc1\tdet2\tloc2\tdet3\tloc3"
//               "\talpha1\tdelta1\taz1\talpha2\tdelta2\taz2\talpha3\tdelta3\taz3\n";
//        for (int i = 0; i < t.size(); i++)
//        {
//            int tpr = t[i].first;
//            auto it1 = std::find_if(f.begin(), f.end(), [tpr](auto& a){return a.first  == tpr;});
//            auto it2 = std::find_if(s.begin(), s.end(), [tpr](auto& a){return a.first  == tpr;});
//            if (it1 != f.end() && it2 != s.end())
//            {
//                double quat3[4] = {t[i].second[0], t[i].second[1], t[i].second[2], t[i].second[3]};
//                double angles3[3];


//                double quat1[4] = {it1->second[0], it1->second[1], it1->second[2], it1->second[3]};
//                double angles1[3];


//                double quat2[4] = {it2->second[0], it2->second[1], it2->second[2], it2->second[3]};
//                double angles2[3];


//                quatToEkvA(quat3, angles3, AxisType::xAxis);
//                quatToEkvA(quat1, angles1, AxisType::yAxis);
//                double Ox31 = calculateScalarProductAngles(angles3, angles1);
//                quatToEkvA(quat2, angles2, AxisType::yAxis);
//                double Ox32 = calculateScalarProductAngles(angles3, angles2);
//                quatToEkvA(quat1, angles1, AxisType::xAxis);
//                double Ox12 = calculateScalarProductAngles(angles2, angles1);



//                //                double Oz31 = calculateAngleAxis(quat3, quat1, AxisType::zAxis);
//                //                double Oz32 = calculateAngleAxis(quat3, quat2, AxisType::zAxis);
//                //                double Oz12 = calculateAngleAxis(quat2, quat1, AxisType::zAxis);
//                quatToEkvA(quat3, angles3);
//                quatToEkvA(quat1, angles1);
//                quatToEkvA(quat2, angles2);
//                double Oz31 = calculateScalarProductAngles(angles3, angles1);
//                double Oz32 = calculateScalarProductAngles(angles3, angles2);
//                double Oz12 = calculateScalarProductAngles(angles2, angles1);

//                //                double Ox31 = calculateAngleAxis(quat3, quat1, AxisType::xAxis);
//                //                double Ox32 = calculateAngleAxis(quat3, quat2, AxisType::xAxis);
//                //                double Ox12 = calculateAngleAxis(quat2, quat1, AxisType::xAxis);

//                out << tpr << "\t" << QString::number(Oz31, 'g', 10) << "\t" <<
//                       QString::number(Oz32,'g', 10) << "\t" <<
//                       QString::number(Oz12,'g', 10) << "\t" <<
//                       QString::number(Ox31,'g', 10) << "\t" <<
//                       QString::number(Ox32,'g', 10) << "\t" <<
//                       QString::number(Ox12,'g', 10) << "\t" <<
//                       it1->second[4] << "\t" << it1->second[5] << "\t" <<
//                       it2->second[4] << "\t" << it2->second[5] << "\t" <<
//                       t[i].second[4] << "\t" << t[i].second[5] << "\t" <<
//                       QString::number(angles1[0],'g', 10) << "\t" << QString::number(angles1[1],'g', 10) << "\t" << QString::number(angles1[2],'g', 10) << "\t" <<
//                                                                                                                                                            QString::number(angles2[0],'g', 10) << "\t" << QString::number(angles2[1],'g', 10) << "\t" << QString::number(angles2[2],'g', 10) << "\t" <<
//                                                                                                                                                                                                                                                                                                 QString::number(angles3[0],'g', 10) << "\t" << QString::number(angles3[1],'g', 10) << "\t" << QString::number(angles3[2],'g', 10) << "\n";
//                //QString::number(angleOzi31,'g', 10) << "\t" << QString::number(angleOzi32,'g', 10) << "\t" << QString::number(angleOzi12,'g', 10) << "\n" ;
//            }
//        }
//    }
