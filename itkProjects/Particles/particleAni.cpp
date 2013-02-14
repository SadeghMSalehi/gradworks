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
#include "piParticleCore.h"



std::vector<pi::Particle> g_Trace;
int g_MaxId = 0;
int g_MaxTimeSteps = 0;
int g_CurrentFrame = 0;
pi::Particle g_BoundingBox;


vtkIntArray* g_IdArray;

AniWindow::AniWindow(QWidget* parent) {
    ui.setupUi(this);
    ui.toolBar->addAction(ui.action_Open);
    ui.toolBar->addAction(ui.actionBackward);
    ui.toolBar->addAction(ui.actionForward);

    m_Renderer = vtkRenderer::New();
    ui.vtkwidget->GetRenderWindow()->AddRenderer(m_Renderer);

    m_Points = vtkPoints::New();
    m_Grid = vtkUnstructuredGrid::New();
    m_Grid->SetPoints(m_Points);

    g_IdArray = vtkIntArray::New();
    g_IdArray->SetNumberOfComponents(1);
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

        m_Renderer->AddViewProp(actor);

    }

    m_Renderer->AddViewProp(actor);
    m_Renderer->ResetCamera();
    ui.vtkwidget->GetRenderWindow()->Render();
}

void AniWindow::UpdateParticles() {
    m_Renderer->RemoveAllViewProps();

    if (g_MaxTimeSteps <= g_CurrentFrame) {
        g_CurrentFrame = g_MaxTimeSteps - 1;
        return;
    } else if (0 > g_CurrentFrame) {
        g_CurrentFrame = 0;
        return;
    }

    for (int i = 0; i < g_MaxId; i++) {
        m_Points->SetPoint(i, g_Trace[g_CurrentFrame*g_MaxId+i].x);
    }

    CreateParticles();

    ui.statusbar->showMessage(QString("%1 Points and %2 Times; Time = %3" ).arg(g_MaxId).arg(g_MaxTimeSteps).arg(g_Trace[g_CurrentFrame*g_MaxId].t));
}

void AniWindow::on_action_Open_triggered() {
    g_Trace.clear();
    m_Points->Reset();
    g_IdArray->Reset();


    FILE* file = fopen("/tmpfs/trace.txt", "r");
    QTextStream in(file);
    QString line;
    for (int cnt = 0; !in.atEnd(); cnt++) {
        pi::Particle p;
        in >> p.t;
        if (in.atEnd()) {
            break;
        }
        in >> p.idx;
        if (p.idx >= g_MaxId) {
            g_MaxId = p.idx + 1;
        }
        g_IdArray->InsertNextValue(p.idx);
        for4(k) in >> p.x[k];
        for4(k) in >> p.y[k];
        for4(k) in >> p.v[k];
        for4(k) in >> p.f[k];
        in >> p.density;
        in >> p.pressure;
        if (cnt == 0) {
            g_BoundingBox = p;
            // x for maximum and y for minimum
            forcopy(p.x, p.y);
        } else {
            formin(g_BoundingBox.x, p.x, g_BoundingBox.x);
            formax(g_BoundingBox.y, p.x, g_BoundingBox.y);
        }
        g_Trace.push_back(p);
    }
    g_MaxTimeSteps = g_Trace.size() / g_MaxId;

    for (int i = 0; i < g_MaxId; i++) {
        m_Points->InsertNextPoint(g_Trace[i].x);
    }

    ui.statusbar->showMessage(QString("%1 Points and %2 Times" ).arg(g_MaxId).arg(g_MaxTimeSteps));

    CreateParticles();
}

void AniWindow::on_actionForward_triggered() {
    if (g_Trace.size() == 0) {
        return;
    }

    g_CurrentFrame ++;
    UpdateParticles();
}


void AniWindow::on_actionBackward_triggered() {
    if (g_Trace.size() == 0) {
        return;
    }

    g_CurrentFrame --;
    UpdateParticles();
}


void AniWindow::on_actionFirst_triggered() {
    if (g_Trace.size() == 0) {
        return;
    }

    g_CurrentFrame = 0;
    UpdateParticles();
}


void AniWindow::on_actionLast_triggered() {
    if (g_Trace.size() == 0) {
        return;
    }

    g_CurrentFrame = g_Trace.size() / g_MaxId - 1;
    UpdateParticles();
}

void AniWindow::on_glyphRadius_valueChanged(double r) {
    UpdateParticles();
}