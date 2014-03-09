//
//  particleAni.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/14/13.
//
//

#include "particleAni.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSphereSource.h"
#include "vtkOutlineSource.h"
#include "vtkGlyph3D.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkPointData.h"
#include "vtkScalarsToColors.h"
#include "vtkLookupTable.h"
#include "vtkActor.h"
#include "vtkIntArray.h"
#include "QTextStream"
#include "QtDebug"
#include "QFile"
#include "QFileDialog"
#include "piParticleCore.h"
#include "piParticleSystemSolver.h"
#include "piParticleTrace.h"
#include "piImageDef.h"
#include "piParticleBSpline.h"
#include "piVTK.h"
#include "piParticleTools.h"


enum { SHOW_TRACE, SHOW_SYSTEM };
static int g_showType = SHOW_TRACE;
static int g_MaxId = 0;
static int g_MaxTimeSteps = 0;
static int g_CurrentFrame = 0;
static int g_CurrentSubject = 0;
static pi::Particle g_BoundingBox;

static pi::ParticleSystemSolver g_Solver;
static pi::ParticleTrace g_Trace;

static vtkIntArray* g_IdArray;
static pivtk::PropScene g_PropScene;

AniWindow::AniWindow(QWidget* parent) {
    ui.setupUi(this);
    ui.toolBar->addAction(ui.actionReset_Camera);
    ui.toolBar->addAction(ui.actionBackward);
    ui.toolBar->addAction(ui.actionForward);
    ui.toolBar->addAction(ui.actionLabel_Registration);

    m_Renderer = vtkRenderer::New();
    ui.vtkwidget->GetRenderWindow()->AddRenderer(m_Renderer);

    m_Points = vtkPoints::New();
    m_Grid = vtkUnstructuredGrid::New();
    m_Grid->SetPoints(m_Points);

    g_IdArray = vtkIntArray::New();
    g_IdArray->SetNumberOfComponents(1);

    if (ui.actionImage_Particle_View->isChecked()) {
        ui.imageDock->show();
    } else {
        ui.imageDock->hide();
    }

    g_PropScene.SetRenderer(m_Renderer);
}

AniWindow::~AniWindow() {
    m_Grid->Delete();
    m_Points->Delete();
}


void AniWindow::CreateParticles() {
    m_Grid->GetPointData()->SetScalars(g_IdArray);

    vtkSphereSource* source = vtkSphereSource::New();
    source->SetRadius(ui.glyphRadius->value());
    source->Update();

    vtkGlyph3D* glyph = vtkGlyph3D::New();
    glyph->SetSource(source->GetOutput());
    glyph->SetInput(m_Grid);
    glyph->SetColorModeToColorByScalar();
    glyph->SetScaleModeToDataScalingOff();
    glyph->Update();

    vtkPolyData* particles = glyph->GetOutput();

    vtkLookupTable* lut = vtkLookupTable::New();
    lut->SetHueRange(0.667, 0.0);
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
    mapper->SetInput(particles);
    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(m_Grid->GetScalarRange());
    mapper->ScalarVisibilityOn();

    vtkActor* actor = vtkActor::New();
    actor->SetMapper(mapper);

    vtkOutlineSource* outline = vtkOutlineSource::New();
    {
        outline->SetBounds(g_BoundingBox.x[0], g_BoundingBox.y[0],
                            g_BoundingBox.x[1], g_BoundingBox.y[1],
                            g_BoundingBox.x[2], g_BoundingBox.y[2]);

        outline->Update();

        vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
        mapper->SetInput(outline->GetOutput());

        vtkActor* actor = vtkActor::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(0.3, 0.3, 0.3);
        actor->GetProperty()->SetOpacity(0.5);

        g_PropScene.AddProp("outline", actor);

    }

    g_PropScene.AddProp("particles", actor);
    ui.vtkwidget->GetRenderWindow()->Render();
}

void AniWindow::OpenTrace(const char* file, bool useFullTrace) {
    g_showType = SHOW_TRACE;
    
    ui.timeSteps->show();
    g_Trace.Clear();
    ui.subjects->clear();
   
    
    ifstream in(file);
    g_Trace.isFullTrace = useFullTrace;
    
    if (in.is_open()) {
        g_Trace.Read(in);
        int n = g_Trace.system.size();
        if (n == 0) {
            return;
        }
        ui.statusbar->showMessage(QString("%1 Points and %2 Times and %3 Subjects" ).arg(g_Trace.system[0].maxIdx + 1).arg(g_Trace.system[0].timeSeries.size()).arg(g_Trace.system.size()));
        
        for (int i = 0; i < n; i++) {
            ui.subjects->addItem(QString("subj-%1").arg(i));
        }
    }
    ui.timeSteps->setMaximum(g_Trace.system[0].timeSeries.size()-1);
    ui.timeStepSpin->setMaximum(g_Trace.system[0].timeSeries.size()-1);
    ui.timeSteps->setEnabled(true);
}

void AniWindow::on_actionOpen_Trace_triggered() {
    QString filename = QFileDialog::getOpenFileName(this, "Open Trace File", ".", "*.txt");
    if (filename.isNull()) {
        return;
    }
    OpenTrace(filename.toUtf8().data(), false);
}

void AniWindow::on_actionOpen_FullTrace_triggered() {
    QString filename = QFileDialog::getOpenFileName(this, "Open Trace File", ".", "*.txt");
    if (filename.isNull()) {
        return;
    }
    OpenTrace(filename.toUtf8().data(), true);
}

void AniWindow::on_actionOpen_System_triggered() {
    g_showType = SHOW_SYSTEM;

    g_MaxTimeSteps = 0;

    m_Points->Reset();
    g_IdArray->Reset();
    ui.subjects->clear();

    QString filename = QFileDialog::getOpenFileName(this, "Open Trace File", ".", "*.txt");
    if (filename.isNull()) {
        return;
    }
    
    if (!g_Solver.LoadConfig(filename.toUtf8().data())) {
        return;
    }

    ui.subjects->addItem("Initial");
    for (int i = 0; i < g_Solver.m_System.GetNumberOfSubjects(); i++) {
        ui.subjects->addItem(QString::fromStdString(g_Solver.m_System[i].m_Name));
    }
    ui.timeSteps->setEnabled(false);
}

void AniWindow::on_subjects_currentIndexChanged(int n) {
    m_Points->Reset();
    g_IdArray->Reset();

    if (g_showType == SHOW_SYSTEM) {
        g_CurrentSubject = n;
        pi::ParticleSubject* subjPtr = &(g_Solver.m_System.GetInitialSubject());
        if (n > 0) {
            subjPtr = &(g_Solver.m_System[n-1]);
            g_CurrentSubject = n - 1;
        }
        pi::ParticleSubject& subj = (*subjPtr);
        for (int i = 0; i < subj.GetNumberOfPoints(); i++) {
            pi::Particle& p = subj[i];
            if (p.idx >= g_MaxId) {
                g_MaxId = p.idx + 1;
            }
            g_IdArray->InsertNextValue(p.idx);
            if (i == 0) {
                g_BoundingBox = p;
                // x for maximum and y for minimum
                forcopy(p.x, p.y);
            } else {
                formin(g_BoundingBox.x, p.x, g_BoundingBox.x);
                formax(g_BoundingBox.y, p.x, g_BoundingBox.y);
            }
            m_Points->InsertNextPoint(p.x);
        }
        g_MaxTimeSteps = 1;
        ui.statusbar->showMessage(QString("%1 Points and %2 Times" ).arg(g_MaxId).arg(g_MaxTimeSteps));

        CreateParticles();

    } else if (g_showType == SHOW_TRACE) {
        g_CurrentSubject = n;
        pi::ParticleSetSeries& snapshot = g_Trace.system[g_CurrentSubject];
        g_MaxTimeSteps = snapshot.timeSeries.size();
        ShowTraceParticles();
    }
}

void AniWindow::ShowTraceParticles() {
    const int showingSubject = ui.subjects->currentIndex();
    pi::ParticleSetSeries& snapshot = g_Trace.system[showingSubject];

    if (g_CurrentFrame < g_MaxTimeSteps) {
        m_Points->Reset();
        g_IdArray->Reset();
        g_BoundingBox = snapshot.boundingBox;
        pi::ParticleVector& data = snapshot.timeSeries[g_CurrentFrame];
        const int npoints = data.size();
        for (int i = 0; i < npoints; i++) {
            g_IdArray->InsertNextValue(i);
            m_Points->InsertNextPoint(data[i].x);
        }
    }
    CreateParticles();
}

void AniWindow::on_actionForward_triggered() {
    if (g_showType != SHOW_TRACE) {
        return;
    }
    if (++g_CurrentFrame >= g_MaxTimeSteps) {
        g_CurrentFrame = g_MaxTimeSteps - 1;
    }
    ShowTraceParticles();
}


void AniWindow::on_actionBackward_triggered() {
    if (g_showType != SHOW_TRACE) {
        return;
    }
    if (--g_CurrentFrame < 0) {
        g_CurrentFrame = 0;
    }
    ShowTraceParticles();
}


void AniWindow::on_actionFirst_triggered() {
    if (g_showType != SHOW_TRACE) {
        return;
    }
    g_CurrentFrame = 0;
    ShowTraceParticles();
}


void AniWindow::on_actionLast_triggered() {
    if (g_showType != SHOW_TRACE) {
        return;
    }
    g_CurrentFrame = g_MaxTimeSteps - 1;
    ShowTraceParticles();
}

void AniWindow::on_actionLabel_Registration_triggered() {
    using namespace pi;

    QString filename = QFileDialog::getOpenFileName(this, "Open Image File", ".", "*.nrrd");
    if (filename.isNull()) {
        return;
    }
    
    itkcmds::itkImageIO<LabelImage> io;
    LabelImage::Pointer sourceImage = io.ReadImageT(filename.toUtf8().data());
    ParticleBSpline particleTransform;
    particleTransform.SetReferenceImage(sourceImage);
    if (g_showType == SHOW_TRACE) {
        if (g_Trace.system.size() > 1) {
            ParticleVector& srcVector = g_Trace.system[0].timeSeries[g_CurrentFrame];
            ParticleVector& dstVector = g_Trace.system[1].timeSeries[g_CurrentFrame];
            particleTransform.EstimateTransform<ParticleXCaster, LabelImage, ParticleVector>(srcVector, dstVector, srcVector.size(), sourceImage);
        } else {
            cout << "Not enough subject for registration" << endl;
            return;
        }
    } else {
        if (g_Solver.m_System.GetNumberOfSubjects() > 1) {
            particleTransform.EstimateTransform(g_Solver.m_System[0], g_Solver.m_System[1]);
        } else {
            cout << "Not enough subject for registration" << endl;
            return;
        }
    }
    LabelImage::Pointer outputImage = particleTransform.WarpLabel(sourceImage);

    QString saveFilename = QFileDialog::getSaveFileName(this, "Save Image File", ".", "*.nrrd");
    if (saveFilename.isNull()) {
        return;
    }
    io.WriteImageT(saveFilename.toUtf8().data(), outputImage);
}


void AniWindow::on_actionMark_At_Image_triggered() {
    QString filename = QFileDialog::getOpenFileName(this, "Open Image File", ".", "*.nrrd");
    if (filename.isNull()) {
        return;
    }

    using namespace pi;
    itkcmds::itkImageIO<LabelImage> io;
    LabelImage::Pointer ref = io.ReadImageT(filename.toUtf8().data());
    LabelImage::Pointer canvas = io.NewImageT(ref);
    ParticleVector& data = g_Trace.system[g_CurrentSubject].timeSeries[g_CurrentFrame];
    MarkAtImage<ParticleVector>(data, data.size(), canvas, 1);

    QString save = QFileDialog::getSaveFileName(this, "Save Image File", ".", "*.nrrd");
    if (save.isNull()) {
        return;
    }
    io.WriteImageT(save.toUtf8().data(), canvas);
}


void AniWindow::on_actionCreate_Pathline_triggered() {
    int n = ui.subjects->currentIndex();
    if (g_showType == SHOW_TRACE) {
        //pivtk::PolyDataPointer data = pivtk::ConstructPathLines(g_Trace.system[n]);
        //g_PropScene.AddPolyData("pathline", data);
        //ui.vtkwidget->GetRenderWindow()->Render();
    }
}

void AniWindow::on_actionReset_Camera_triggered() {
    m_Renderer->ResetCamera();
    ui.vtkwidget->GetRenderWindow()->Render();
}

void AniWindow::on_actionImage_Particle_View_triggered(bool checked) {
    if (checked) {
        ui.imageDock->show();
    } else {
        ui.imageDock->hide();
    }
}

void AniWindow::on_timeStepSpin_valueChanged(int value) {
    ui.timeSteps->setValue(ui.timeStepSpin->value());
}


void AniWindow::on_timeSteps_valueChanged(int value) {
    if (g_showType == SHOW_TRACE) {
        g_CurrentFrame = ui.timeSteps->value();
        ui.timeStepSpin->setValue(g_CurrentFrame);
        ShowTraceParticles();
    }
}

void AniWindow::on_glyphRadius_valueChanged(double r) {
    on_subjects_currentIndexChanged(ui.subjects->currentIndex());
}