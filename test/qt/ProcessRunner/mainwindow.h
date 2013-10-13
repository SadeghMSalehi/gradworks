#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>
#include <string>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    std::string findProgram(std::string name);

public slots:
    void runScript();
    void stopScript();
    void readScriptOutput();
    void scriptStateChanged(QProcess::ProcessState state);

    void addRow();
    void deleteCurrentRow();
    void editCell(int row, int col);

private:
    Ui::MainWindow *ui;
    QProcess _script;
};

#endif // MAINWINDOW_H
