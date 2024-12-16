### RK_solver.c 
Programa para solucionar la ecuación de van der Pol mediante el método explicito Runge-Kuta de cuarto orden,
el programa utiliza los 5 parámetros de ejecución, y registra la solución para 3 seg en archivos de texto
nombrados y_sol.txt, para la solución con respecto a y, y el archivo phase_data.txt de las 2 variables y, u.

*Para observar las graficas resultantes es necesario tener instalado GNUPLOT.

Los argumentos de ejecución son:
*todos son obligatorios
1. V1 (flotante o entero)
2. V2 (flotante o entero de signo opuesto a v1 para obtener oscilador relajado)
3. mu (flotante o entero no negativo)
4. e (flotante o entero)
5. d (flotante o entero)


### vdpApp.cpp
Programa que implementa el programa descrito en RK_solver.c haciendo uso de una interfaz grafica con wxWidgets en C++,
es necesario tener todas las dependencias wxWidgets así como de GNUPLOT.

Pre-requisitos:
 Comandos de Linux utilizados para la configuración de wxWidgets desde un
 sub-sistema en windows (WSL) en su distibución de Ubuntu:

 1. Instalación de compilador (si no se encuentra con el comando $g++ --version):

    $sudo apt update

    $sudo apt install build-essential

 2. Instalación de librerias GTK para graficos de wxWidgets:

    $sudo apt install libgtk-3-dev

    $sudo apt install libgl1-mesa-dev

    $sudo apt install libglu1-mesa-dev

    $sudo apt install libwebkit2gtk-4.0-dev

 3. Instalación de wxWidgets desde el repositorio basado en la documentación disponible en https://github.com/wxWidgets/wxWidgets/blob/master/docs/gtk/install.md
   Los comandos suponen una distribución basada en Debian tal como lo es Ubuntu:

    $git clone --recurse-submodules https://github.com/wxWidgets/wxWidgets.git

    $cd wxWidgets/

    $mkdir buildgtk

    $cd buildgtk/

    $../configure --with-gtk

    (El numero sobre el comando siguiente coincide con los nucleos en el procesador que ejecutará el programa,
    esto puede consultarse con el comando "$nproc")

    $make -j16 

    
    //Para verificar que el proceso de instalación ha sido correcto hasta este punto, se puede ejecutar alguno de los codigos de prueba incluidos 

    //en los archivos instalados de wxWidgets, para el ejemplo más sencillo esto es:

        $cd samples/minimal

        $make

        $./minimal


### Ejecución del proyecto:
Con la carpeta del proyecto localizada en una ruta X del sistema, se puede ejecutar a partir del archivo "Makefile" como:

    $cd X/proyectoVDP/build

    $make

    $./vdpApp
    
//En caso de que se muestren advertencias estilo 

    Gdk-Message: Unable to load sb_v_double_arrow from the cursor theme

//se puede ejecutar el siguiente comando para obtener el estilo de graficos que solicita el programa

    $sudo apt-get -y install adwaita-icon-theme-full
