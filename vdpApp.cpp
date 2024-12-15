/* Programa basado en wxWidgets dedicado a la solución
    de la ecuación de van der Pol utilizando el metodo 
    de Runge-Kuta de cuarto orden.

    Métodos Numéricos 
    M. en C. de Computacion
    Edgar Iván Aguilera Hernández
    12 / 12 /2024
*/

#include "wx/wx.h"
#include <wx/filedlg.h>
#include "wx/sizer.h"
#include "wx/rawbmp.h"
#include "wx/splitter.h"
#include "wx/spinctrl.h"
#include <iostream>

/*van der pol equation plot*/
class vdpSolver{
    private:
        //equation system parameters
        double range[2] = {0,3};        //solution range 
        static const int n = 500;       //number of steps 
        static const int m = 2;         //number of equations/variables
        double init_c[2] = {2.0,0.0};   //initial conditions y(0), v(0)

        double v[2];                    //depolarization - fractal speed 
        double a;                       //damping factor 
        double d;
        double e; 
        double solution[n+1][m+1];
        
        
        //Van Der Pol equation as 2 first grade equations with variables t,y,u
        double vdpOscillator(double *var, int equation){

            double solution;
            if(equation == 0){
                // y'= u
                solution = var[1];
            }else if(equation == 1){
                //u' = k(1-y^2)u - ay
                solution = a * (1 - var[0]*var[0])*var[1] - d*var[0];
            }
            return solution;
        }

        //Van Der Pol modified equation as 2 first grade equations with variables t,y,u
        double vdpOscillatorModified(double *var, int equation){
            double result;

            if(equation == 0){
                // y'= v
                result = var[1];
            }else if(equation == 1){
                //v' = -k(y - v_1)(y  - v_2)u - [y(y + d)(y + e)] / de
                result = -a * (var[0]-v[0])*(var[0]-v[1] )*var[1] - 
                (var[0]*(var[0] + d)*(var[0] + e)) / d*e;
            }

            return result;
        }

        //4th order Runge-Kuta method for m ordinary differential equations
        void solveRungeKuta(){
            
            //set of 4th grade approximations, m equations
            double k[4][m];
            //set size of steps
            double h = (range[1] - range[0]) /n;
            //set starting value 
            double t = range[0];
            solution[0][0] = t;
            //set initial solution as the initial condition
            double w[m];
            double w_aux[m];
            copyVector(init_c,w,m);
            for(int i = 0; i< m; i++){
                solution[0][i+1] = w[i];
            }
            

            //iterate solution approximation over the n steps of [a,b]
            for(int i = 0; i < n; i++){
                //start series of averages to reach 4th-order
                //k-1
                for(int j = 0; j < m; j++){
                    k[0][j] = h*vdpOscillatorModified(w,j);
                }
                //add previous estimations to the variables (excluding indepedent variable)
                for(int l = 0; l < m; l++){
                    w_aux[l] = w[l] + (0.5)*k[0][l]; 
                } 
                //k-2
                for(int j = 0; j < m; j++){
                    k[1][j] = h*vdpOscillatorModified(w_aux,j);
                }
                //add previous estimations to the variables (excluding indepedent variable)
                for(int l = 0; l < m; l++){
                    w_aux[l] = w[l] + (0.5)*k[1][l]; 
                } 
                //k-3
                for(int j = 0; j < m; j++){
                    k[2][j] = h*vdpOscillatorModified(w_aux,j);
                }
                //add previous estimations to the variables (excluding indepedent variable)
                for(int l = 0; l < m; l++){
                    w_aux[l] = w[l] + k[2][l]; 
                } 
                //k-4
                for(int j = 0; j < m; j++){
                    k[3][j] = h*vdpOscillatorModified(w_aux,j);
                }

                //update solution of w with k-th estimations
                for(int j = 0; j < m; j++){
                    w[j] += (k[0][j] + 2*k[1][j] + 2*k[2][j] + k[3][j])/6.0;
                }
                
                //update t with next step
                t += h; 

                //store solution set of variables for t value
                solution[i+1][0] = t;
                for(int j = 0; j< m; j++){
                    solution[i+1][j+1] = w[j];
                }
            }

        }

        //copy double vector 1(1xn) into 2(1xn) 
        void copyVector(double *vector_1, double *vector_2, int n){
            for(int i=0; i<n; i++){
                vector_2[i] = vector_1[i];
            }
        }

        void plotGNU(){
            //iniciar aplicacion GNUPLOT
            FILE *fPlotter = popen("gnuplot -persist", "w");
            
            //Configurar plot
            fprintf(fPlotter,"set terminal png size 1000,650\n");
            fprintf(fPlotter,"set title \"Potencial de acción: v_1: %.1f, v_2: %.1f, a: %1.f, d: %1.f, e: %1.f\" \n",
                    v[0],v[1],a,d,e);
            fprintf(fPlotter,"set xlabel \" tiempo(s) \" \n");
            fprintf(fPlotter,"set xrange [%f:%f] \n",range[0],range[1]);
            fprintf(fPlotter,"set ylabel \" U \" \n");
            fprintf(fPlotter,"set yrange [-3:3] \n");
            fprintf(fPlotter,"set grid ytics mytics \n");
            fprintf(fPlotter,"set grid \n");

            //crear grafica
            fprintf(fPlotter, "plot \"/home/edgarbanzo/Programs/NM/proyectoVDP/resources/y_sol.txt\" using 1:2 smooth csplines lw 3 lt 6 lc rgb \"red\" title \"Nodo Sinoatrial \" \n");
            fprintf(fPlotter,"set output \"/home/edgarbanzo/Programs/NM/proyectoVDP/resources/vdp_curve%d.png\"\n",n);
            fprintf(fPlotter,"replot\n");

            //phase diagram
            fprintf(fPlotter,"set title \"Plano de fase\" \n");
            fprintf(fPlotter,"set xlabel \" y(mv) \" \n");
            fprintf(fPlotter,"set xrange [-4:4] \n");
            fprintf(fPlotter,"set ylabel \" u \" \n");
            fprintf(fPlotter,"set yrange [-15:15] \n");

            fprintf(fPlotter, "plot \"/home/edgarbanzo/Programs/NM/proyectoVDP/resources/phase_data.txt\" with lines lw 3 lt 2 title \"Orbita de Sol.\"\n");
            fprintf(fPlotter,"set terminal png size 1000,650\n");
            fprintf(fPlotter,"set output \"/home/edgarbanzo/Programs/NM/proyectoVDP/resources/phase_diagram%d.png\"\n",n);
            fprintf(fPlotter,"replot\n");
            
            fclose(fPlotter);
        }

    public:
        //constructor set values to rest state parameters
        vdpSolver(){
            //solutions matrix (n+1 solutions, m variables + independent variable)
            for(int i = 0; i < n+1; i++){
                for(int j = 0; j < m+1; j++){
                    solution[i][j] = 0.0;
                }
            }
            v[0] = 0.8;
            v[1] = -0.8;
            a = 15.0; 
            d = 6.0;
            e = 8.0;

        }

        void setParameters(double v1_new,double v2_new, double a_new, double d_new, double e_new){
            v[0] = v1_new;
            v[1] = -v2_new;
            a = a_new;
            d = d_new;
            e = e_new;
        }

        void setv1(double v1_new){v[0] = v1_new;}
        void setv2(double v2_new){v[1] = -v2_new;}
        void seta(double a_new){a = a_new;}
        void setd(double d_new){d = d_new;}
        void sete(double e_new){e = e_new;}

        //solve van der pole equation returning path to the plotted solution
        wxString getPlot(int equation){
            //solve ordinal differential equation
            solveRungeKuta();

            //write solution into file to plot
            //write data points to plot
            FILE *fData_y = fopen("/home/edgarbanzo/Programs/NM/proyectoVDP/resources/y_sol.txt","w");
            FILE *fData_p = fopen("/home/edgarbanzo/Programs/NM/proyectoVDP/resources/phase_data.txt","w");
            for(int i = 0; i < n+1; i++){
                fprintf(fData_y,"%0.4f %0.4f\n",solution[i][0],solution[i][1]);
                fprintf(fData_p,"%0.4f %0.4f\n",solution[i][1],solution[i][2]);
            }
            fclose(fData_y);
            fclose(fData_p);

            //plot interpolated data
            plotGNU();

            wxString plotPath = wxString::Format(wxT("/home/edgarbanzo/Programs/NM/proyectoVDP/resources/vdp_curve%d.png"),n);

            return plotPath;
        }
};


/*Image panel handler*/
class wxImagePanel : public wxPanel
{
    wxImage image;
    wxBitmap m_bitmap;
    wxBitmap resized;
    int w, h;
    
public:
    wxImagePanel(wxSplitterWindow *parent, wxString file, wxBitmapType format);
    wxImagePanel(wxSplitterWindow *parent);
    void setImage(wxString file, wxBitmapType format);
    void paintEvent(wxPaintEvent & evt);
    void paintNow();
    void OnSize(wxSizeEvent& event);
    void render(wxDC& dc);
    //static event handling
    DECLARE_EVENT_TABLE();
};


BEGIN_EVENT_TABLE(wxImagePanel, wxPanel)
    EVT_PAINT(wxImagePanel::paintEvent)// catch paint events
    EVT_SIZE(wxImagePanel::OnSize)//Size event
END_EVENT_TABLE()

/*constructor with default local image (test purposes)*/
wxImagePanel::wxImagePanel(wxSplitterWindow *parent, wxString file, wxBitmapType format) :
wxPanel(parent){
    
    m_bitmap.LoadFile(file, format);
    image = wxImage(m_bitmap.ConvertToImage());
    w = image.GetWidth();
    h = image.GetHeight();
}

/*starting program constructor*/
wxImagePanel::wxImagePanel(wxSplitterWindow *parent):wxPanel(parent){
    //default black image of size 100 x 100
    image = wxImage(100,100,true);
    w = image.GetWidth();
    h = image.GetHeight();
}

/*set new image loaded from path*/
void wxImagePanel::setImage(wxString file, wxBitmapType format){
    
    m_bitmap.LoadFile(file, format);
    image = wxImage(m_bitmap.ConvertToImage());
    w = image.GetWidth();
    h = image.GetHeight();
}

/*Refresh the image whith any change or event associated (triggered manually by calling Refresh()/Update())*/
void wxImagePanel::paintEvent(wxPaintEvent & evt){
    wxPaintDC dc(this);
    render(dc);
}

void wxImagePanel::paintNow(){
    wxClientDC dc(this);
    render(dc);
}

void wxImagePanel::render(wxDC&  dc){
    int neww, newh;
    dc.GetSize( &neww, &newh );
    
    if( neww != w || newh != h )
    {
        resized = wxBitmap( image.Scale( neww, newh /*, wxIMAGE_QUALITY_HIGH*/ ) );
        w = neww;
        h = newh;
        dc.DrawBitmap( resized, 0, 0, false );
    }else{
        dc.DrawBitmap( resized, 0, 0, false );
    }
}

/*tell the panel to draw itself again (when the user resizes the image panel)*/
void wxImagePanel::OnSize(wxSizeEvent& event){
    Refresh();
    event.Skip();
}

 
class MyFrame : public wxFrame{
    public:
        MyFrame(wxBoxSizer *sizer);

    private:
        vdpSolver vdpEquation;
        wxImagePanel *drawPanel;    //instance of panel where graph is displayed
        wxImagePanel *gifPanel;     //instance of panel where heart is displayed
        wxPanel *optionPanel;       //instance of panel where options are displayed
        wxComboBox *statesList;     //pointer to instance of combobox
        wxSpinCtrlDouble *v_1;      //pointer to instance of "v1" spin control
        wxSpinCtrlDouble *v_2;      //pointer to instance of "v2" spin control
        wxSpinCtrlDouble *d;        //pointer to instance of "d" spin control
        wxSpinCtrlDouble *e;        //pointer to instance of "e" spin control
        wxSpinCtrlDouble *a;        //pointer to instance of "a" spin control
        wxButton *apply;            //pointer to instance of "aplicar" button
        wxTimer *timer;             //pointer to instance of timer to change gif frame
        int frameCount = 0;             //frame id to display

        void resetFrame();

        //static event handling
        void OnComboBoxSelection(wxCommandEvent& event);
        void OnButtonApplyClick(wxCommandEvent& event);
        void OnV1SpinChange(wxCommandEvent& event);
        void OnV2SpinChange(wxCommandEvent& event);
        void OnaSpinChange(wxCommandEvent& event);
        void OndSpinChange(wxCommandEvent& event);
        void OneSpinChange(wxCommandEvent& event);
        void ChangeFrame(wxTimerEvent& event);
        DECLARE_EVENT_TABLE();
        
};
 

/*Element's label assignation*/
enum{
    MPLOT = 2,
    BUTTON1 = 3,
    COMBOBOX = 5,
    SPINCTRL1 = 6,
    SPINCTRL2 = 7,
    SPINCTRL3 = 8,
    SPINCTRL4 = 9,
    SPINCTRL5 = 11,
};

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_COMBOBOX(COMBOBOX,MyFrame::OnComboBoxSelection)
    EVT_BUTTON(BUTTON1,MyFrame::OnButtonApplyClick)
END_EVENT_TABLE()

//////////////////////////////////////////////////////////////////window elements initialization


/*Elemnts initialization associated with the main window*/
MyFrame::MyFrame(wxBoxSizer *sizer)
    : wxFrame(nullptr, wxID_ANY, "Proyecto van der Pol I",wxPoint(1200,1200), wxSize(1200,700)){
    
    //add splitters for better UX
    wxSplitterWindow *splitter = new wxSplitterWindow(this,wxID_ANY, wxDefaultPosition,wxDefaultSize,
                                wxSP_BORDER | wxSP_LIVE_UPDATE);
    
    wxSplitterWindow *rightsplitter = new wxSplitterWindow(splitter,wxID_ANY, wxDefaultPosition,wxDefaultSize,
                                wxSP_BORDER | wxSP_LIVE_UPDATE);
    
    drawPanel = new wxImagePanel(splitter);
    gifPanel = new wxImagePanel(rightsplitter);          
    optionPanel = new wxPanel(rightsplitter);

    sizer = new wxBoxSizer(wxHORIZONTAL);
    splitter->SetSizer(sizer);
    sizer->Add(drawPanel, 1, wxEXPAND);
    wxBoxSizer *gifSizer = new wxBoxSizer(wxVERTICAL);
    rightsplitter->SetSizer(gifSizer);
    gifSizer->Add(gifPanel,1,wxEXPAND | wxALL, 5);

    splitter->SetMinimumPaneSize(950);
    splitter->SplitVertically(drawPanel,rightsplitter);
    rightsplitter->SetMinimumPaneSize(340);
    rightsplitter->SplitHorizontally(optionPanel,gifPanel);
    rightsplitter->SetSashPosition(0);
    rightsplitter->SetSashGravity(1);

    //initialize gif timer
    timer = new wxTimer();

    //buttons initialization
    //apply = new wxButton(optionPanel,BUTTON1,_T("Aplicar"),wxPoint(160,45));
    //apply->SetBackgroundColour(wxColour(117, 240, 230));
    wxString choices[] = {_T("Reposo"),
                        _T("Taquicardia"),
                        _T("Bradicardia")};
    statesList = new wxComboBox(optionPanel,COMBOBOX,"Reposo",wxPoint(10,45), wxSize(140,30),
                        3, choices);
    
    v_1 = new wxSpinCtrlDouble(optionPanel,SPINCTRL1,"0.8",wxPoint(35,120),wxSize(129,35),
                    wxSP_ARROW_KEYS,0.0,10.0,1.0,0.2);
    v_2 = new wxSpinCtrlDouble(optionPanel,SPINCTRL2,"0.8",wxPoint(35,160),wxSize(129,35),
                    wxSP_ARROW_KEYS,0.0,10.0,1.0,0.2);
    a = new wxSpinCtrlDouble(optionPanel,SPINCTRL3,"15.0",wxPoint(32,210),wxSize(129,35),
                    wxSP_ARROW_KEYS,0.0,20.0,1.0,1.0); 
    d = new wxSpinCtrlDouble(optionPanel,SPINCTRL4,"6.0",wxPoint(32,250),wxSize(129,35),
                    wxSP_ARROW_KEYS,0.0,20.0,1.0,1.0);
    e = new wxSpinCtrlDouble(optionPanel,SPINCTRL5,"8.0",wxPoint(32,290),wxSize(129,35),
                    wxSP_ARROW_KEYS,0.0,20.0,1.0,1.0);
    
    
    //static text initialization
    new wxStaticText(optionPanel,wxID_ANY,"Estados:",wxPoint(10,20),wxDefaultSize);
    new wxStaticText(optionPanel,wxID_ANY,"Parametros:",wxPoint(10,100),wxDefaultSize);
    new wxStaticText(optionPanel,wxID_ANY,"V1:",wxPoint(10,129),wxDefaultSize);
    new wxStaticText(optionPanel,wxID_ANY,"V2:",wxPoint(10,169),wxDefaultSize);
    new wxStaticText(optionPanel,wxID_ANY,"a:",wxPoint(10,219),wxDefaultSize);
    new wxStaticText(optionPanel,wxID_ANY,"d:",wxPoint(10,259),wxDefaultSize);
    new wxStaticText(optionPanel,wxID_ANY,"e:",wxPoint(10,299),wxDefaultSize);

    //Status Message at the bottom of the window
    CreateStatusBar();
    SetStatusText("Proyecto Metodos numericos I, V_1.0");
    //van der Pol equation solver initialization
    vdpEquation = vdpSolver();

    //dynamic events binding
    Bind(wxEVT_SPINCTRLDOUBLE,&MyFrame::OnV1SpinChange,this,SPINCTRL1);
    Bind(wxEVT_SPINCTRLDOUBLE,&MyFrame::OnV2SpinChange,this,SPINCTRL2);
    Bind(wxEVT_SPINCTRLDOUBLE,&MyFrame::OnaSpinChange,this,SPINCTRL3);
    Bind(wxEVT_SPINCTRLDOUBLE,&MyFrame::OndSpinChange,this,SPINCTRL4);
    Bind(wxEVT_SPINCTRLDOUBLE,&MyFrame::OneSpinChange,this,SPINCTRL5);
    timer->Bind(wxEVT_TIMER, &MyFrame::ChangeFrame, this);
    timer->Start(80);
}


/*reset values on spin controls to their default states*/
void MyFrame::resetFrame(){

}
 

/*Selection made on type of processing list*/
void MyFrame::OnComboBoxSelection(wxCommandEvent& event){

    
}

/*Changes made on V1 spin controls to change equation parameter*/
void MyFrame::OnV1SpinChange(wxCommandEvent& event){
    //set new equation parameters 
    vdpEquation.setv1(v_1->GetValue());
    //solve equation
    drawPanel->setImage(vdpEquation.getPlot(1),wxBITMAP_TYPE_PNM);
    drawPanel->Refresh();
}

/*Changes made on V2 spin controls to change equation parameter*/
void MyFrame::OnV2SpinChange(wxCommandEvent& event){
    //set new equation parameters 
    vdpEquation.setv2(v_2->GetValue());
    //solve equation
    drawPanel->setImage(vdpEquation.getPlot(1),wxBITMAP_TYPE_PNM);
    drawPanel->Refresh();
}

/*Changes made on a spin controls to change equation parameter*/
void MyFrame::OnaSpinChange(wxCommandEvent& event){
    //set new equation parameters 
    vdpEquation.seta(a->GetValue());
    //solve equation
    drawPanel->setImage(vdpEquation.getPlot(1),wxBITMAP_TYPE_PNM);
    drawPanel->Refresh();
}

/*Changes made on d spin controls to change equation parameter*/
void MyFrame::OndSpinChange(wxCommandEvent& event){
    //set new equation parameters 
    vdpEquation.setd(d->GetValue());
    //solve equation
    drawPanel->setImage(vdpEquation.getPlot(1),wxBITMAP_TYPE_PNM);
    drawPanel->Refresh();
}

/*Changes made on e spin controls to change equation parameter*/
void MyFrame::OneSpinChange(wxCommandEvent& event){
    //set new equation parameters 
    vdpEquation.sete(e->GetValue());
    //solve equation
    drawPanel->setImage(vdpEquation.getPlot(1),wxBITMAP_TYPE_PNM);
    drawPanel->Refresh();
}

/*Click on "Aplicar" button */
void MyFrame::OnButtonApplyClick(wxCommandEvent& event){

    //set new equation parameters 
    vdpEquation.setParameters(v_1->GetValue(),v_2->GetValue(),a->GetValue(),d->GetValue(),e->GetValue());
    //solve equation
    drawPanel->setImage(vdpEquation.getPlot(1),wxBITMAP_TYPE_PNM);
    drawPanel->Refresh();

}

/*Change gif frame on timer event*/
void MyFrame::ChangeFrame(wxTimerEvent& event){
    //reset gif
    if(frameCount > 17)
        frameCount = 0;
    
    wxString framePath = wxString::Format(wxT("/home/edgarbanzo/Programs/NM/proyectoVDP/resources/heartGif/frame_%d.png")
                                        ,frameCount);
    gifPanel->setImage(framePath,wxBITMAP_TYPE_PNG);
    gifPanel->Refresh();
    frameCount++;
}


////////////////////////////////////////////////////////////////////////Application Initialization (wxwidgets main)
class MyApp : public wxApp{

public:
    bool OnInit() override;//virtual function to initialize application
};
 
wxIMPLEMENT_APP(MyApp);
 
bool MyApp::OnInit(){
    wxInitAllImageHandlers();
    wxBoxSizer *sizer;
    MyFrame *frame = new MyFrame(sizer);
    frame->Show(true);
    return true;
}