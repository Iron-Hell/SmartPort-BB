function simulacion_piezoelectrico
    % =====================================================================
    % SIMULACIÓN DE RECOLECTOR DE ENERGÍA POR IMPACTO (PÉNDULO)
    % MODELO: CABECEO (PITCH) CON ROTACIÓN RESPECTO AL CENTRO DE FLOTACIÓN
    % =====================================================================
    
    clear; clc; close all;

    % --- PARÁMETROS GEOMÉTRICOS Y FÍSICOS ---
    w = 0.10;          % Ancho del recinto (10 cm)
    h = 0.15;          % Altura del recinto (15 cm)
    L = h / 2;         % Longitud del péndulo (7.5 cm)
    m = 0.05;          % Masa del péndulo (50 gramos)
    g = 9.81;          % Aceleración de la gravedad
    
    % --- CINEMÁTICA REAL DE LA BOYA ---
    % Distancia vertical desde el centro de flotación (CF) hasta el eje de anclaje del péndulo.
    R_piv = 0.20;      % 20 cm por debajo del anclaje
    
    % --- PROPIEDADES DEL SISTEMA ---
    zeta = 0.02;       % Coeficiente de amortiguamiento del péndulo
    r = 0.65;          % Coeficiente de restitución mecánica en el impacto piezoeléctrico
    eta = 0.20;        % Eficiencia de conversión electromecánica del material (20%)
    omega_n = sqrt(g / L); % Frecuencia natural angular del péndulo
    
    % --- EXCITACIÓN EXTERNA (MOVIMIENTO DE LAS OLAS) ---
    f_ext = 2.0;       % Frecuencia del cabeceo de la ola en Hz (2 Hz)
    omega_ext = 2 * pi * f_ext;
    theta_0 = 5 * (pi / 180); % Amplitud máxima de cabeceo (5 grados convertidos a radianes)
    
    % --- LÍMITES DE IMPACTO (PAREDES) ---
    % Ángulo límite relativo en radianes donde se encuentran los parches piezoeléctricos
    phi_lim = asin((w / 2) / L); 
    
    % --- PARÁMETROS DE TIEMPO DE LA SIMULACIÓN ---
    t_sim = 100;       % Ventana de tiempo representativo a simular (segundos)
    t_start = 0;
    y0 = [0; 0];       % Condiciones iniciales: [ángulo_relativo (rad), velocidad_angular (rad/s)]
    
    % --- VARIABLES DE ALMACENAMIENTO ---
    t_todo = [];
    y_todo = [];
    energia_total_cosechada = 0;
    num_colisiones = 0;
    
    % Configuración de las opciones del solucionador con detección de eventos (choques)
    options = odeset('Events', @(t, y) evento_colision(t, y, phi_lim), 'RelTol', 1e-6, 'AbsTol', 1e-8);
    
    % --- BUCLE PRINCIPAL DE INTEGRACIÓN (MÉTODO PARTICIONADO POR IMPACTOS) ---
    while t_start < t_sim
        % Resolver numéricamente hasta que ocurra una colisión o termine el tiempo
        [t, y, te, ye, ie] = ode45(@(t, y) ecuaciones_pendulo(t, y, omega_n, zeta, theta_0, omega_ext, R_piv, L), ...
                                   [t_start, t_sim], y0, options);
        
        % Acumular los datos de la trayectoria actual
        t_todo = [t_todo; t]; 
        y_todo = [y_todo; y]; 
        
        % Comprobar si la integración se detuvo por un evento de colisión
        if ~isempty(te)
            num_colisiones = num_colisiones + 1;
            t_start = te(end); % Avanzar el tiempo inicial al instante del choque
            
            % Velocidad angular justo antes del impacto (rad/s)
            vel_impacto = ye(end, 2); 
            
            % Energía cinética transferida / perdida en el choque
            E_mecanica_perdida = 0.5 * m * (L * vel_impacto)^2 * (1 - r^2);
            
            % Energía eléctrica aprovechada según la eficiencia del piezoeléctrico
            energia_total_cosechada = energia_total_cosechada + (eta * E_mecanica_perdida);
            
            % Aplicar la condición de rebote instantáneo (velocidad invertida y atenuada)
            y0 = [ye(end, 1); -r * vel_impacto];
        else
            % Se completó el tiempo de simulación sin más colisiones laterales
            break; 
        end
    end
    
    % --- EXTRAPOLACIÓN A 24 HORAS ---
    segundos_24h = 24 * 60 * 60;          % Segundos totales en un día (86400 s)
    factor_escala = segundos_24h / t_sim;  % Factor para escalar de t_sim a 24h
    
    colisiones_24h = num_colisiones * factor_escala;
    energia_24h = energia_total_cosechada * factor_escala;
    potencia_promedio = energia_total_cosechada / t_sim; % Potencia media continua (Watts)
    
    % --- REPORTE DE RESULTADOS EN CONSOLA ---
    fprintf('=====================================================\n');
    fprintf('         RESULTADOS DEL MODELO PREDICTIVO            \n');
    fprintf('=====================================================\n');
    fprintf('Tiempo de análisis simulado       : %d segundos\n', t_sim);
    fprintf('Colisiones registradas en ventana : %d impactos\n', num_colisiones);
    fprintf('-----------------------------------------------------\n');
    fprintf('Colisiones estimadas en 24 horas  : %.0f impactos\n', colisiones_24h);
    fprintf('Energía cosechada en 24 horas     : %.4f Joules\n', energia_24h);
    fprintf('Potencia promedio del dispositivo : %.4f mWatts\n', potencia_promedio * 1000);
    fprintf('=====================================================\n');
    
    % --- GRÁFICA DE CONTROL DE LA DINÁMICA ---
    figure('Color', 'w');
    plot(t_todo, y_todo(:,1) * (180/pi), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ángulo del péndulo'); 
    hold on;
    yline(phi_lim * (180/pi), 'r--', 'Piezoeléctrico x=w', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');
    yline(-phi_lim * (180/pi), 'r--', 'Piezoeléctrico x=0', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');
    
    xlabel('Tiempo (s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Ángulo relativo \phi (grados)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Respuesta Temporal del Péndulo con Efecto de Brazo de Palanca Flotante', 'FontSize', 12);
    grid on;
    xlim([0 min(t_sim, 20)]); % Mostrar solo los primeros 20 segundos para apreciar los choques con claridad
    legend('Location', 'best');
end

% =====================================================================
% SUBSISTEMA: ECUACIONES DIFERENCIALES DEL PÉNDULO (MARCO NO INERCIAL)
% =====================================================================
function dydt = ecuaciones_pendulo(t, y, omega_n, zeta, theta_0, omega_ext, R_piv, L)
    phi = y(1);   % Posición angular relativa
    dphi = y(2);  % Velocidad angular relativa
    
    % Aceleración angular de cabeceo (derivada segunda de theta)
    % ddtheta = -theta_0 * (omega_ext^2) * sin(omega_ext * t);
    % Debido al sistema de referencia rotativo no inercial, el signo cambia al pasar como fuerza ficticia.
    % Además se incorpora el efecto de acoplamiento por el brazo de palanca al pivote real (R_piv).
    
    factor_acoplamiento = 1 + (R_piv / L);
    excitacion_inercial = factor_acoplamiento * theta_0 * (omega_ext^2) * sin(omega_ext * t);
    
    % Ecuación de movimiento diferencial ordinaria:
    ddphi = -2 * zeta * omega_n * dphi - (omega_n^2) * phi + excitacion_inercial;
    
    dydt = [dphi; ddphi];
end

% =====================================================================
% SUBSISTEMA: DETECTOR DE EVENTOS DE COLISIÓN MECÁNICA
% =====================================================================
function [value, isterminal, direction] = evento_colision(~, y, phi_lim)
    % Evaluamos la cercanía al límite físico en valor absoluto.
    % Cuando "value" llega a cero, significa que tocó la pared izquierda o derecha.
    value = phi_lim - abs(y(1)); 
    
    isterminal = 1;  % Detener inmediatamente la integración para aplicar condiciones de rebote
    direction = 0;   % Detectar impacto sin importar si va de izquierda a derecha o viceversa
end