#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <mlstask.h>
#include <QScopedPointer>
#include <QMessageBox>
#include <QSettings>
#include <QDebug>
#include <algorithm>
#include <customplot/cxyplotter.h>
#include <QMessageBox>

struct PointDiff
{
    QCPGraph* graph;
    QCPItemTracer* tracer;
};

using namespace std;
using namespace BOKZMath;
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_chooseRawFilesPushButton_clicked();

    void on_saveToolButton_clicked();

    void on_removeToolButton_clicked();

    void on_startToolButton_clicked();

private:

    bool checkChoose();

    void saveSettings();

    void loadSettings();

    void printResults(Results& res);

    void printErrors(ResultErrors& err);

    void saveErrors(ResultErrors& err);

    void printAngles (const QString& before, const QString& after);

    void saveResults(Results& res);

    QBitArray setFlags();

    Results setFirstApprox();


    Ui::MainWindow *ui;
    QScopedPointer <CXYPlotter> plotter;
    QVector <PointDiff> diffPlotInfo;
    QString fileName;
    MLSTask task;
    QSettings* settings = nullptr;
    bool editStarted = false;
    QStringList editingList;
    qint32 selectedIndex = 0;
};

#endif // MAINWINDOW_H
