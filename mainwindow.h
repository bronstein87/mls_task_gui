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
using namespace std;
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_chooseRawFilesPushButton_clicked();

private:
    void joinMeasureFiles();
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
    QString fileName;
    MLSTask task;
    QSettings* settings = nullptr;
};

#endif // MAINWINDOW_H
