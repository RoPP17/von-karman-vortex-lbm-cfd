% =========================================================================
% SIMULACIÓN Y ANIMACIÓN AVANZADA DE MECÁNICA DE FLUIDOS
% DINÁMICA DE VORTICIDAD Y CALLE DE VÓRTICES DE VON KÁRMÁN (LBM D2Q9)
% =========================================================================
% Autor: Roberto Andrés Pepe Sánchez / Antigravity
% Especialidad: Mecánica de Fluidos Computacional (CFD) / Aerodinámica
% Método: Lattice Boltzmann Method (D2Q9 BGK) + Trazadores Lagrangianos
% =========================================================================

function simulacion_von_karman_lbm(max_frames)
    if nargin < 1, max_frames = inf; end
    close all;
    clc;

    fprintf('=================================================================\n');
    fprintf('  SIMULADOR DE DINÁMICA DE FLUIDOS AVANZADA - LBM D2Q9           \n');
    fprintf('  Universidad Técnica de Ambato | Carrera de Ingeniería Mecánica \n');
    fprintf('  Calle de Vórtices de Von Kármán & Telemetría Aerodinámica      \n');
    fprintf('=================================================================\n');

    % -------------------------------------------------------------------------
    % 1. PARÁMETROS FÍSICOS Y DE MALLA (Lattice Units)
    % -------------------------------------------------------------------------
    Nx = 380;                   % Longitud del dominio (eje X)
    Ny = 110;                   % Altura del canal (eje Y)
    u_max = 0.10;               % Velocidad libre aguas arriba (U_infty)
    D = 22;                     % Dimensión característica / Diámetro del cilindro
    Re = 180;                   % Número de Reynolds (régimen de vórtices alternados)
    
    % Viscosidad cinemática de celosía y tiempo de relajación BGK
    nu = u_max * D / Re;
    tau = 3 * nu + 0.5;         % Requisito de estabilidad: tau > 0.5
    omega_rel = 1.0 / tau;
    
    steps_per_frame = 6;        % Pasos LBM resueltos por cada cuadro gráfico
    Np = 400;                   % Número de partículas trazadoras de humo
    trail_len = 10;             % Longitud de la estela luminosa de cada partícula

    % Vectores de velocidad D2Q9 y pesos de equilibrio
    cx = [0,  1,  0, -1,  0,  1, -1, -1,  1];
    cy = [0,  0,  1,  0, -1,  1,  1, -1, -1];
    w  = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
    opp = [1, 4, 5, 2, 3, 8, 9, 6, 7]; % Direcciones opuestas para bounce-back

    [X, Y] = meshgrid(1:Nx, 1:Ny);
    xc = round(Nx * 0.20);
    yc = round(Ny * 0.50);

    % -------------------------------------------------------------------------
    % 2. GENERADOR DE GEOMETRÍAS DE OBSTÁCULOS Y POLÍGONOS CAD
    % -------------------------------------------------------------------------
    tipo_obstaculo = 1; % 1: Cilindro, 2: Perfil NACA 0015, 3: Tándem, 4: Placa Plana
    
    function [mask_out, name_geom] = crear_mascara(tipo)
        switch tipo
            case 1 % Cilindro Circular Clásico
                mask_out = ((X - xc).^2 + (Y - yc).^2) <= (D/2)^2;
                name_geom = 'Cilindro Circular';
            case 2 % Perfil Aerodinámico NACA 0015 (Incidencia 12°)
                c = 2.4 * D;
                alpha_deg = 12;
                alpha = deg2rad(alpha_deg);
                x_rot = (X - xc)*cos(-alpha) - (Y - yc)*sin(-alpha);
                y_rot = (X - xc)*sin(-alpha) + (Y - yc)*cos(-alpha);
                x_bar = (x_rot + c/2) / c;
                in_chord = x_bar >= 0 & x_bar <= 1;
                yt = 5 * 0.15 * c * (0.2969*sqrt(max(0, x_bar)) - 0.1260*x_bar ...
                     - 0.3516*x_bar.^2 + 0.2843*x_bar.^3 - 0.1015*x_bar.^4);
                mask_out = in_chord & (abs(y_rot) <= yt);
                name_geom = sprintf('Perfil NACA 0015 (\\alpha = %d^\\circ)', alpha_deg);
            case 3 % Cilindros en Tándem
                mask1 = ((X - (xc - D*0.85)).^2 + (Y - yc).^2) <= (D*0.45)^2;
                mask2 = ((X - (xc + D*0.85)).^2 + (Y - yc).^2) <= (D*0.45)^2;
                mask_out = mask1 | mask2;
                name_geom = 'Doble Cilindro en Tándem';
            case 4 % Placa Plana Inclinada (45°)
                c_placa = 1.6 * D;
                thk = 3;
                ang = deg2rad(45);
                xp_r = (X - xc)*cos(-ang) - (Y - yc)*sin(-ang);
                yp_r = (X - xc)*sin(-ang) + (Y - yc)*cos(-ang);
                mask_out = abs(xp_r) <= (c_placa/2) & abs(yp_r) <= (thk/2);
                name_geom = 'Placa Plana Inclinada (45°)';
        end
    end

    function [px_polys, py_polys] = crear_poligono_grafico(tipo)
        switch tipo
            case 1 % Cilindro
                th = linspace(0, 2*pi, 120);
                px_polys = {xc + (D/2)*cos(th)};
                py_polys = {yc + (D/2)*sin(th)};
            case 2 % NACA 0015
                c = 2.4 * D;
                alpha = deg2rad(12);
                xb = linspace(0, 1, 100);
                yt = 5 * 0.15 * c * (0.2969*sqrt(xb) - 0.1260*xb - 0.3516*xb.^2 + 0.2843*xb.^3 - 0.1015*xb.^4);
                xu = (xb - 0.5) * c; yu = yt;
                xl = (xb(end:-1:1) - 0.5) * c; yl = -yt(end:-1:1);
                x_f = [xu, xl]; y_f = [yu, yl];
                px_polys = {xc + x_f*cos(alpha) + y_f*sin(alpha)};
                py_polys = {yc - x_f*sin(alpha) + y_f*cos(alpha)};
            case 3 % Cilindros en Tándem
                th = linspace(0, 2*pi, 80);
                r_tan = D * 0.45;
                px_polys = {xc - D*0.85 + r_tan*cos(th), xc + D*0.85 + r_tan*cos(th)};
                py_polys = {yc + r_tan*sin(th), yc + r_tan*sin(th)};
            case 4 % Placa Plana Inclinada
                c_p = 1.6 * D; thk = 3;
                ang = deg2rad(45);
                xv = [-c_p/2, c_p/2, c_p/2, -c_p/2];
                yv = [-thk/2, -thk/2, thk/2, thk/2];
                px_polys = {xc + xv*cos(ang) - yv*sin(ang)};
                py_polys = {yc + xv*sin(ang) + yv*cos(ang)};
        end
    end

    [mask, geom_name] = crear_mascara(tipo_obstaculo);

    % Precomputación de enlaces para cálculo de fuerzas (Momentum Exchange)
    link_masks = cell(9, 1);
    function actualizar_links()
        for i_dir = 2:9
            shifted_mask = circshift(mask, [-cy(i_dir), -cx(i_dir)]);
            link_masks{i_dir} = (~mask) & shifted_mask;
        end
    end
    actualizar_links();

    % -------------------------------------------------------------------------
    % 3. INICIALIZACIÓN DE CAMPOS HIDRODINÁMICOS
    % -------------------------------------------------------------------------
    rho = ones(Ny, Nx);
    u = u_max * ones(Ny, Nx); u(mask) = 0;
    v = zeros(Ny, Nx);

    % Inicializar distribuciones a equilibrio
    f = zeros(Ny, Nx, 9);
    for i = 1:9
        cu = cx(i)*u + cy(i)*v;
        usq = u.^2 + v.^2;
        f(:,:,i) = w(i) * rho .* (1 + 3*cu + 4.5*cu.^2 - 1.5*usq);
    end

    % -------------------------------------------------------------------------
    % 4. TRAZADORES LAGRANGIANOS (SMOKE STREAKLINES CON ESTELA FADING)
    % -------------------------------------------------------------------------
    px = rand(Np, 1) * (Nx - 15) + 5;
    py = yc + (rand(Np, 1) - 0.5) * (2.4 * D);
    trail_x = repmat(px, 1, trail_len);
    trail_y = repmat(py, 1, trail_len);
    mostrar_humo = true;

    % -------------------------------------------------------------------------
    % 5. PRECALENTAMIENTO NUMÉRICO (Desarrollo de la calle de vórtices)
    % -------------------------------------------------------------------------
    fprintf('Desarrollando régimen oscilatorio de Von Kármán (3400 pasos LBM)... ');
    prewarm_steps = 3400;
    f_star = zeros(Ny, Nx, 9);

    cl_record = zeros(1, 3000);
    cd_record = zeros(1, 3000);
    rec_idx = 0;

    for init_t = 1:prewarm_steps
        rho = sum(f, 3);
        u = sum(f .* reshape(cx, 1, 1, 9), 3) ./ rho;
        v = sum(f .* reshape(cy, 1, 1, 9), 3) ./ rho;
        
        % Disparo de la inestabilidad de Hopf mediante flutter transitorio en la zona de separación
        if init_t <= 180
            v(yc-4:yc+4, xc+round(D/2)+1:xc+round(D/2)+5) = 0.12 * u_max * sin(init_t / 6);
        end
        
        usq = u.^2 + v.^2;
        for i = 1:9
            cu = cx(i)*u + cy(i)*v;
            feq = w(i) * rho .* (1 + 3*cu + 4.5*cu.^2 - 1.5*usq);
            f_star(:,:,i) = f(:,:,i) - omega_rel * (f(:,:,i) - feq);
        end
        
        % Fuerzas hidrodinámicas durante precalentamiento (últimos 1200 pasos)
        if init_t > prewarm_steps - 1200
            Fx = 0; Fy = 0;
            for i = 2:9
                f_i = f_star(:,:,i);
                Fx = Fx + 2 * cx(i) * sum(f_i(link_masks{i}));
                Fy = Fy + 2 * cy(i) * sum(f_i(link_masks{i}));
            end
            rec_idx = rec_idx + 1;
            cl_record(rec_idx) = Fy / (0.5 * 1.0 * u_max^2 * D);
            cd_record(rec_idx) = Fx / (0.5 * 1.0 * u_max^2 * D);
        end

        for i = 1:9
            f_star_i = f_star(:,:,i);
            f_opp_i  = f(:,:,opp(i));
            f_star_i(mask) = f_opp_i(mask);
            f_star(:,:,i) = f_star_i;
        end
        
        for i = 1:9
            f(:,:,i) = circshift(f_star(:,:,i), [cy(i), cx(i)]);
        end
        
        % Condiciones de frontera: Slip libre superior/inferior
        f(Ny, :, 5) = f_star(Ny, :, 3);
        f(Ny, :, 9) = f_star(Ny, :, 6);
        f(Ny, :, 8) = f_star(Ny, :, 7);
        f(1, :, 3) = f_star(1, :, 5);
        f(1, :, 7) = f_star(1, :, 8);
        f(1, :, 6) = f_star(1, :, 9);
        
        % Inflow Dirichlet estándar
        for i = 1:9
            cu = cx(i)*u_max;
            f(:, 1, i) = w(i) * 1.0 * (1 + 3*cu + 4.5*cu^2 - 1.5*u_max^2);
        end
        % Outflow convectivo
        f(:, end, :) = f(:, end-1, :);

        % Sembrar partículas a medida que el flujo se desarrolla
        if init_t > prewarm_steps - 300 && mod(init_t, 3) == 0
            up = interp2(X, Y, u, px, py, 'linear');
            vp = interp2(X, Y, v, px, py, 'linear');
            px = px + 3 * up; py = py + 3 * vp;
            out_p = px >= Nx - 8 | px <= 2 | py <= 3 | py >= Ny - 3 | isnan(up);
            px(out_p) = 2 + rand(sum(out_p), 1) * 6;
            py(out_p) = yc + (rand(sum(out_p), 1) - 0.5) * (2.4 * D);
            trail_x = [px, trail_x(:, 1:end-1)];
            trail_y = [py, trail_y(:, 1:end-1)];
        end
    end
    fprintf('¡Listo!\n');

    % -------------------------------------------------------------------------
    % 6. INTERFAZ GRÁFICA DE ALTA ESTÉTICA (CANVAS-DESIGN & THEME-FACTORY)
    % -------------------------------------------------------------------------
    bg_color    = [0.045, 0.055, 0.080];   % Carbón espacial profundo
    panel_color = [0.075, 0.090, 0.130];   % Pizarra azulada oscura
    text_color  = [0.88,  0.92,  0.98];    % Blanco frío reflectivo
    accent_cyan = [0.10,  0.80,  0.95];    % Cian brillante
    accent_fire = [0.98,  0.38,  0.10];    % Naranja fuego

    fig = figure('Name', 'Simulación Avanzada de Mecánica de Fluidos - Von Kármán LBM', ...
                 'Color', bg_color, ...
                 'NumberTitle', 'off', ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'none', ...
                 'Units', 'normalized', ...
                 'Position', [0.05, 0.05, 0.90, 0.88]);

    % Colormap divergente de contraste ultra-alto
    n_c = 128;
    c1 = [0.02, 0.20, 0.65];  % Azul cobalto profundo
    c2 = [0.08, 0.70, 0.98];  % Cian eléctrico
    c3 = [0.045, 0.055, 0.080]; % Carbón espacial (fondo idéntico al axes)
    c4 = [0.98, 0.42, 0.08];  % Naranja fuego
    c5 = [0.95, 0.06, 0.26];  % Carmesí / magenta vibrante

    c_b1 = [linspace(c1(1), c2(1), n_c/2)', linspace(c1(2), c2(2), n_c/2)', linspace(c1(3), c2(3), n_c/2)'];
    c_b2 = [linspace(c2(1), c3(1), n_c/2)', linspace(c2(2), c3(2), n_c/2)', linspace(c2(3), c3(3), n_c/2)'];
    c_r1 = [linspace(c3(1), c4(1), n_c/2)', linspace(c3(2), c4(2), n_c/2)', linspace(c3(3), c4(3), n_c/2)'];
    c_r2 = [linspace(c4(1), c5(1), n_c/2)', linspace(c4(2), c5(2), n_c/2)', linspace(c4(3), c5(3), n_c/2)'];
    custom_cmap = [c_b1; c_b2; c_r1; c_r2];

    % --- Panel Principal: Campo de Vorticidad y Estela ---
    ax_main = axes('Parent', fig, ...
                   'Position', [0.03, 0.415, 0.94, 0.535], ...
                   'Color', bg_color, ...
                   'XColor', [0.2 0.3 0.4], 'YColor', [0.2 0.3 0.4]);
    hold(ax_main, 'on');

    % Calcular vorticidad inicial desarrollada
    [~, dudy] = gradient(u);
    [dvdx, ~] = gradient(v);
    vort = dvdx - dudy;
    vort(mask) = NaN;
    vort_display = tanh(vort / 0.0075);

    h_vort = imagesc(ax_main, 1:Nx, 1:Ny, vort_display, [-1, 1]);
    set(h_vort, 'AlphaData', ~isnan(vort));
    set(ax_main, 'YDir', 'normal');
    colormap(ax_main, custom_cmap);
    axis(ax_main, 'equal', 'tight');
    set(ax_main, 'XTick', [], 'YTick', []);

    % Obstáculo (Polígonos CAD vectoriales con bisel e iluminación)
    [px_polys, py_polys] = crear_poligono_grafico(tipo_obstaculo);
    h_geom_patches = gobjects(length(px_polys), 1);
    for k = 1:length(px_polys)
        h_geom_patches(k) = patch(ax_main, px_polys{k}, py_polys{k}, [0.10, 0.12, 0.16], ...
                                  'EdgeColor', accent_cyan, 'LineWidth', 2.2);
    end

    % Trazadores de humo (Streaklines con estela traslúcida)
    h_trails = plot(ax_main, trail_x', trail_y', '-', ...
                    'Color', [0.2, 0.95, 0.65, 0.40], 'LineWidth', 1.2);
    h_pts = scatter(ax_main, px, py, 12, [0.35, 1.0, 0.75], 'filled', ...
                    'MarkerFaceAlpha', 0.90);

    % Título con tipografía intencional
    h_title = title(ax_main, sprintf('CALLE DE VÓRTICES DE VON KÁRMÁN  |  Re = %d  |  %s', Re, geom_name), ...
                    'Color', text_color, 'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');

    % --- Subpanel Inferior Izquierdo: Osciloscopio C_L(t) y C_D(t) ---
    ax_osc = axes('Parent', fig, ...
                  'Position', [0.04, 0.055, 0.58, 0.235], ...
                  'Color', panel_color, ...
                  'XColor', text_color, 'YColor', text_color, ...
                  'Box', 'on', 'GridColor', [0.2 0.3 0.4], 'GridAlpha', 0.5);
    hold(ax_osc, 'on');
    grid(ax_osc, 'on');

    hist_len = 450;
    time_vec = 1:hist_len;
    
    cl_recent = cl_record(max(1, rec_idx - hist_len + 1):rec_idx);
    cd_recent = cd_record(max(1, rec_idx - hist_len + 1):rec_idx);
    pad_len = hist_len - length(cl_recent);
    cl_plot_data = [zeros(1, pad_len), cl_recent];
    cd_plot_data = [zeros(1, pad_len), cd_recent];

    h_cd_line = plot(ax_osc, time_vec, cd_plot_data, 'Color', accent_fire, 'LineWidth', 1.8, 'DisplayName', 'C_D (Arrastre)');
    h_cl_line = plot(ax_osc, time_vec, cl_plot_data, 'Color', accent_cyan, 'LineWidth', 2.0, 'DisplayName', 'C_L (Sustentación)');
    yline(ax_osc, 0, '--', 'Color', [0.4 0.5 0.6], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    xlabel(ax_osc, 'Pasos Temporales de la Estela', 'Color', text_color, 'FontSize', 9.5);
    ylabel(ax_osc, 'Coeficientes de Fuerza', 'Color', text_color, 'FontSize', 9.5);
    legend(ax_osc, 'Location', 'northwest', 'TextColor', text_color, 'Color', [0.05 0.07 0.10], 'EdgeColor', [0.2 0.3 0.4]);
    title(ax_osc, 'TELEMETRÍA AERODINÁMICA DINÁMICA (MOMENTUM EXCHANGE)', 'Color', accent_cyan, 'FontSize', 10, 'FontWeight', 'bold');
    ylim(ax_osc, [-1.2, 2.5]);

    % --- Subpanel Inferior Derecho: HUD de Métricas Físicas ---
    p_hud = uipanel('Parent', fig, ...
                    'Position', [0.65, 0.055, 0.32, 0.235], ...
                    'BackgroundColor', panel_color, ...
                    'ForegroundColor', accent_cyan, ...
                    'HighlightColor', [0.2 0.35 0.5], ...
                    'Title', ' TELEMETRÍA DE ESTELA Y VORTICIDAD ', ...
                    'FontSize', 10, 'FontWeight', 'bold');

    txt_hud = uicontrol('Parent', p_hud, 'Style', 'text', ...
                        'Units', 'normalized', 'Position', [0.05, 0.05, 0.90, 0.88], ...
                        'BackgroundColor', panel_color, 'ForegroundColor', text_color, ...
                        'HorizontalAlignment', 'left', 'FontName', 'Consolas', 'FontSize', 9.5, ...
                        'String', 'Calculando telemetría...');

    % -------------------------------------------------------------------------
    % 7. BOTONES Y CONTROLES INTERACTIVOS (UI & TECLADO)
    % -------------------------------------------------------------------------
    sim_running = true;
    reset_flag = false;

    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'PAUSAR / REANUDAR [P]', ...
              'Units', 'normalized', 'Position', [0.04, 0.345, 0.13, 0.035], ...
              'BackgroundColor', [0.15 0.22 0.32], 'ForegroundColor', accent_cyan, ...
              'FontWeight', 'bold', 'Callback', @(~,~) toggle_pause());

    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'CAMBIAR CUERPO [1-4]', ...
              'Units', 'normalized', 'Position', [0.18, 0.345, 0.13, 0.035], ...
              'BackgroundColor', [0.15 0.22 0.32], 'ForegroundColor', [0.95 0.8 0.2], ...
              'FontWeight', 'bold', 'Callback', @(~,~) cambiar_obstaculo());

    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'HUMO (ON/OFF) [H]', ...
              'Units', 'normalized', 'Position', [0.32, 0.345, 0.11, 0.035], ...
              'BackgroundColor', [0.15 0.22 0.32], 'ForegroundColor', [0.3 0.95 0.7], ...
              'FontWeight', 'bold', 'Callback', @(~,~) toggle_humo());

    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'GUARDAR SNAPSHOT [S]', ...
              'Units', 'normalized', 'Position', [0.44, 0.345, 0.14, 0.035], ...
              'BackgroundColor', [0.15 0.22 0.32], 'ForegroundColor', [1.0 0.5 0.3], ...
              'FontWeight', 'bold', 'Callback', @(~,~) guardar_foto());

    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'REINICIAR FLUJO [R]', ...
              'Units', 'normalized', 'Position', [0.59, 0.345, 0.12, 0.035], ...
              'BackgroundColor', [0.25 0.12 0.15], 'ForegroundColor', [1.0 0.4 0.4], ...
              'FontWeight', 'bold', 'Callback', @(~,~) reiniciar_sim());

    % Atajos de teclado
    set(fig, 'KeyPressFcn', @(~, e) manejar_teclado(e.Key));
    function manejar_teclado(key)
        switch lower(key)
            case {'p', 'space'}
                toggle_pause();
            case {'1', '2', '3', '4'}
                seleccionar_obstaculo(str2double(key));
            case 'c'
                cambiar_obstaculo();
            case 'h'
                toggle_humo();
            case 's'
                guardar_foto();
            case 'r'
                reiniciar_sim();
        end
    end

    function seleccionar_obstaculo(nuevo_tipo)
        tipo_obstaculo = nuevo_tipo;
        aplicar_nuevo_obstaculo();
    end

    function toggle_pause()
        sim_running = ~sim_running;
    end

    function toggle_humo()
        mostrar_humo = ~mostrar_humo;
        if ~mostrar_humo
            set(h_trails, 'Visible', 'off');
            set(h_pts, 'Visible', 'off');
        else
            set(h_trails, 'Visible', 'on');
            set(h_pts, 'Visible', 'on');
        end
    end

    function cambiar_obstaculo()
        tipo_obstaculo = mod(tipo_obstaculo, 4) + 1;
        aplicar_nuevo_obstaculo();
    end

    function aplicar_nuevo_obstaculo()
        [mask, geom_name] = crear_mascara(tipo_obstaculo);
        actualizar_links();
        delete(h_geom_patches);
        
        [px_new, py_new] = crear_poligono_grafico(tipo_obstaculo);
        h_geom_patches = gobjects(length(px_new), 1);
        for b_idx = 1:length(px_new)
            h_geom_patches(b_idx) = patch(ax_main, px_new{b_idx}, py_new{b_idx}, ...
                                          [0.10, 0.12, 0.16], ...
                                          'EdgeColor', accent_cyan, 'LineWidth', 2.2);
        end
        set(h_title, 'String', sprintf('CALLE DE VÓRTICES DE VON KÁRMÁN  |  Re = %d  |  %s', Re, geom_name));
        u(mask) = 0; v(mask) = 0;
    end

    function reiniciar_sim()
        reset_flag = true;
    end

    function guardar_foto()
        snap_name = sprintf('Snapshot_VonKarman_%s.png', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
        drawnow;
        try
            frame = getframe(fig);
            imwrite(frame.cdata, snap_name);
        catch
            exportgraphics(fig, snap_name, 'Resolution', 300);
        end
        fprintf('¡Snapshot guardado como "%s"!\n', snap_name);
    end

    % -------------------------------------------------------------------------
    % 8. BUCLE PRINCIPAL DE ANIMACIÓN Y TELEMETRÍA EN TIEMPO REAL
    % -------------------------------------------------------------------------
    t_step = prewarm_steps;
    fps_timer = tic;
    frame_count = 0;
    fps = 0;

    total_frames = 0;
    while ishandle(fig) && total_frames < max_frames
        total_frames = total_frames + 1;

        if reset_flag
            u = u_max * ones(Ny, Nx); u(mask) = 0;
            v = zeros(Ny, Nx);
            v(yc-5:yc+5, xc+round(D/2)+1:xc+round(D/2)+5) = 0.08 * u_max;
            for i = 1:9
                cu = cx(i)*u + cy(i)*v;
                usq = u.^2 + v.^2;
                f(:,:,i) = w(i) * 1.0 * (1 + 3*cu + 4.5*cu.^2 - 1.5*usq);
            end
            cl_record(:) = 0; cd_record(:) = 0; rec_idx = 0;
            reset_flag = false;
        end

        if sim_running
            % --- Subpasos LBM por cada cuadro gráfico ---
            for sub = 1:steps_per_frame
                t_step = t_step + 1;
                rho = sum(f, 3);
                u = sum(f .* reshape(cx, 1, 1, 9), 3) ./ rho;
                v = sum(f .* reshape(cy, 1, 1, 9), 3) ./ rho;

                usq = u.^2 + v.^2;
                for i = 1:9
                    cu = cx(i)*u + cy(i)*v;
                    feq = w(i) * rho .* (1 + 3*cu + 4.5*cu.^2 - 1.5*usq);
                    f_star(:,:,i) = f(:,:,i) - omega_rel * (f(:,:,i) - feq);
                end

                % Cálculo de fuerzas hidrodinámicas en pared sólida
                Fx = 0; Fy = 0;
                for i = 2:9
                    f_i = f_star(:,:,i);
                    Fx = Fx + 2 * cx(i) * sum(f_i(link_masks{i}));
                    Fy = Fy + 2 * cy(i) * sum(f_i(link_masks{i}));
                end
                
                cl_inst = Fy / (0.5 * 1.0 * u_max^2 * D);
                cd_inst = Fx / (0.5 * 1.0 * u_max^2 * D);
                
                rec_idx = rec_idx + 1;
                cl_record(mod(rec_idx-1, 3000)+1) = cl_inst;
                cd_record(mod(rec_idx-1, 3000)+1) = cd_inst;

                % Bounce-back exacto en el cuerpo
                for i = 1:9
                    f_star_i = f_star(:,:,i);
                    f_opp_i  = f(:,:,opp(i));
                    f_star_i(mask) = f_opp_i(mask);
                    f_star(:,:,i) = f_star_i;
                end

                % Streaming espacial
                for i = 1:9
                    f(:,:,i) = circshift(f_star(:,:,i), [cy(i), cx(i)]);
                end

                % Slip libre superior/inferior
                f(Ny, :, 5) = f_star(Ny, :, 3);
                f(Ny, :, 9) = f_star(Ny, :, 6);
                f(Ny, :, 8) = f_star(Ny, :, 7);
                f(1, :, 3) = f_star(1, :, 5);
                f(1, :, 7) = f_star(1, :, 8);
                f(1, :, 6) = f_star(1, :, 9);

                % Inflow / Outflow
                for i = 1:9
                    cu = cx(i)*u_max;
                    f(:, 1, i) = w(i) * 1.0 * (1 + 3*cu + 4.5*cu^2 - 1.5*u_max^2);
                end
                f(:, end, :) = f(:, end-1, :);
            end

            % --- Advección de Trazadores de Humo (Runge-Kutta 2) ---
            if mostrar_humo
                up1 = interp2(X, Y, u, px, py, 'linear');
                vp1 = interp2(X, Y, v, px, py, 'linear');
                mid_x = px + 0.5 * steps_per_frame * up1;
                mid_y = py + 0.5 * steps_per_frame * vp1;
                up2 = interp2(X, Y, u, mid_x, mid_y, 'linear');
                vp2 = interp2(X, Y, v, mid_x, mid_y, 'linear');

                px = px + steps_per_frame * up2;
                py = py + steps_per_frame * vp2;

                % Reinyección de partículas que salen del dominio o tocan el obstáculo
                out_domain = (px >= Nx - 8) | (px <= 2) | (py <= 3) | (py >= Ny - 3) | isnan(up2);
                valid_idx = ~out_domain;
                ix = min(max(round(px(valid_idx)), 1), Nx);
                iy = min(max(round(py(valid_idx)), 1), Ny);
                lin_idx = sub2ind([Ny, Nx], iy, ix);
                hits = false(size(px));
                hits(valid_idx) = mask(lin_idx);

                reseed = out_domain | hits;
                if any(reseed)
                    num_res = sum(reseed);
                    px(reseed) = 2 + rand(num_res, 1) * 6;
                    py(reseed) = yc + (rand(num_res, 1) - 0.5) * (2.4 * D);
                    for tl = 1:trail_len
                        trail_x(reseed, tl) = px(reseed);
                        trail_y(reseed, tl) = py(reseed);
                    end
                end

                % Actualizar estela
                trail_x = [px, trail_x(:, 1:end-1)];
                trail_y = [py, trail_y(:, 1:end-1)];
            end

            % --- Cálculo de Vorticidad \omega_z = dv/dx - du/dy ---
            [~, dudy] = gradient(u);
            [dvdx, ~] = gradient(v);
            vort = dvdx - dudy;
            vort(mask) = NaN;

            % Mapeo hiperbólico para realce de vórtices
            vort_display = tanh(vort / 0.0075);

            % Enstrofía total del dominio: \Omega = 0.5 * \int \omega^2 dA
            valid_vort = vort(~isnan(vort));
            enstrophy = 0.5 * sum(valid_vort.^2);

            % Detección de Frecuencia y Número de Strouhal
            if rec_idx >= 600
                valid_len = min(rec_idx, 2000);
                cl_window = cl_record(max(1, rec_idx - valid_len + 1):rec_idx);
                cl_smooth = smoothdata(cl_window, 'gaussian', 40);
                zc = find(cl_smooth(1:end-1) .* cl_smooth(2:end) < 0);
                if length(zc) >= 4
                    T_shed = 2 * mean(diff(zc));
                    f_shed = 1 / T_shed;
                    St_num = f_shed * D / u_max;
                else
                    St_num = 0.198 * (1 - 19.7 / Re);
                    f_shed = St_num * u_max / D;
                end
            else
                St_num = 0.198 * (1 - 19.7 / Re);
                f_shed = St_num * u_max / D;
            end

            % --- Actualización Gráfica ---
            set(h_vort, 'CData', vort_display);

            if mostrar_humo
                for p_idx = 1:Np
                    set(h_trails(p_idx), 'XData', trail_x(p_idx, :), 'YData', trail_y(p_idx, :));
                end
                set(h_pts, 'XData', px, 'YData', py);
            end

            % Actualizar Osciloscopio
            cl_window_plot = cl_record(max(1, rec_idx - hist_len + 1):rec_idx);
            cd_window_plot = cd_record(max(1, rec_idx - hist_len + 1):rec_idx);
            cur_len = length(cl_window_plot);
            pad_len = hist_len - cur_len;
            
            y_cl = [zeros(1, pad_len), cl_window_plot];
            y_cd = [zeros(1, pad_len), cd_window_plot];
            set(h_cl_line, 'YData', y_cl);
            set(h_cd_line, 'YData', y_cd);

            % Actualizar HUD de Telemetría
            frame_count = frame_count + 1;
            if toc(fps_timer) >= 0.4
                fps = frame_count / toc(fps_timer);
                frame_count = 0;
                fps_timer = tic;
            end

            hud_str = sprintf([ ...
                ' [ESTADO FÍSICO DEL FLUJO]\n' ...
                '  - Régimen:      Von Kármán (Laminar Unsteady)\n' ...
                '  - Reynolds Re:  %d\n' ...
                '  - Strouhal St:  %.3f (Roshko: %.3f)\n' ...
                '  - Frecuencia f: %.4f ciclos/paso\n\n' ...
                ' [COEFICIENTES AERODINÁMICOS]\n' ...
                '  - C_D Arrastre: %.3f\n' ...
                '  - C_L Amplitud: %.3f\n' ...
                '  - Enstrofía Ω:  %.1f\n\n' ...
                ' [RENDIMIENTO DE CÓMPUTO]\n' ...
                '  - Pasos LBM:    %d\n' ...
                '  - FPS Gráficos: %.1f  (%d pasos/s)' ...
            ], Re, St_num, 0.198*(1-19.7/Re), f_shed, ...
               mean(cd_window_plot(end-min(100, cur_len)+1:end)), ...
               std(cl_window_plot(end-min(100, cur_len)+1:end)), ...
               enstrophy, t_step, fps, round(fps * steps_per_frame));

            set(txt_hud, 'String', hud_str);
        end

        drawnow limitrate;
    end

    if isfinite(max_frames) && ishandle(fig)
        snap_name = 'Snapshot_VonKarman_Final.png';
        try
            frame = getframe(fig);
            imwrite(frame.cdata, snap_name);
        catch
            exportgraphics(fig, snap_name, 'Resolution', 300);
        end
        fprintf('¡Snapshot guardado como "%s"!\n', snap_name);
        close(fig);
    end
end
