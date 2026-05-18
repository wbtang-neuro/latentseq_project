classdef SweepsSettings

    properties

        % Some handy context variables
        verbose = 0
        numericClass = "single"
        useGpu = false  % if true, cast() will convert all grid vectors to gpuArrays
        angleUnits = "degrees"

        % dataRoot_ = "\\forskning.it.ntnu.no\ntnu\mh-kin\moser\sweeps\richarga\temp\sweeps_data\sharing_v4"
        dataRoot_ = 'E:\Tang_latentseq_NN_2026\sample_data\'
        useLocalFileCache = true
        codeRoot_ = [] % leave empty to use default
        useRelativeDataPaths = false

        lockDatasets = 0 % write-protect original dataset files
    end

    properties (Dependent)
        gve       % grid edge vectors
        gv        % grid center vectors
        gg        % mesh-grids in true shape
        gg1       % vectorized mesh-grids
        gvrange   % grid-vector ranges
        gvspacing
        gverange  % grid-vector edge ranges
        info
    end

    properties (Constant)

        col_id = [0.1, 0.6, 0.1]
        col_hd = [0.2, 0.2, 0.8]
        col_theta = [0.5, 0.5, 0.5]
        col_covmodel = [0.3, 0.5, 1]

        col_left = [1, 0, 0];
        col_right = [0, 0, 1];

        col_pos_true = [0.5, 0.5, 0.5]
        col_pos_sweep = [0, 0, 0]

        col_conj = [.9, 0.1, .6];
        col_conjunctive = [.9, 0.1, .6];
        col_grid = [.1, .1, 1];
        col_nongrid = [0, 0, 0] + 0.3
        col_grid_m1 = [0.85, 0, 0]
        col_grid_m2 = [0.3, 0, 1]
        col_grid_m3 = [0.6, 0.6, 0]
        col_grid_m4 = [0, 0, 0]

        col_celltypes = struct( ...
            "b", [0.2,  0.2,  1.0 ], ...
            "c", [1,    0.5,  0   ], ...
            "n", [0.2,  0.6,  0.2 ], ...
            "p", [0.75, 0,    0.75], ...
            "u", [0.6,  0.6,  0.6 ], ...
            "i", [0.0,  0.0,  0.0 ], ...
            "g", [0.2,  0.2,  1.0 ]);


        col_cyc_even = [0.9961, 0.5255, 0.1216];
        col_cyc_odd = [0.2000, 0.6000, 0.8510];

        colormap = viridis()

        cmaps = struct( ...
            "firing_rate", viridis(), ...
            "probability", hot(256), ...
            "connectivity", seismic(256));

        recs_of_mec = ["25843_1"    "25691_2"    "24365_2"    "26035_1"    "25953_4"    "26034_3"    "25127_1"    "26018_2"    "26820_2"    "26648_1"    "27764_1"    "27765_1"    "28063_4"    "28304_1"    "28258_4"    "28229_2"]
        recs_of_mec_hc = ["26035_1"    "26034_3"    "28063_4"    "28258_4"    "28229_2"    "29502_1"];
        recs_ww = ["25843_1" "25691_1"];
        recs_lt = ["26648_1"    "27765_2"    "28063_1"    "28229_3"    "28304_1"    "29502_3"]
        recs_novel
        recs_sleep


        lmt_pos_lenffR = 5
        lmt_pos_lenffR_gmod = [5.5, 4, 3]
        lmt_chunk_length = 1e4
        lmt_types = ["id+pos"]
        lmt_brain_area_combs = {"mec", "hc", ["mec", "hc"]}
        valid_locations = ["mec", "hc"]

        minSpeed = 0.15; % Used for sweep analysis
        max_contamination_rate = 0.3; % pretty lax

        ratemap_nbins_close = 1

        pos_lims = struct( ...
            "of",            [-0.75, 0.75], ...
            "of_1m",         [-0.5, 0.5], ...
            "ww",            [-0.9, 0.9], ...
            "lt",            [-0.75, 0.75], ...
            "mmaze",         [-0.7, 0.7] );

        pos_nbins = struct( ...
            "coarse",        15, ...
            "fine",          60); % gives 2.5 cm bin width

        pos_extend_factor = 1.6 % multiplier for bin count
        
        gve_ = struct( ...
            "angular", linspace(-pi, pi, 60+1)', ...
            "theta_interp", linspace(0, 1, 60+1)', ...
            "torus", linspace(-pi, pi, 25+1)');


        pos_ndilate = 2

        dt = 10e-3

        tsm = struct("pos", 1, "pos_slow", 15, "id", 1, "speed", 100); % sigma in samples

        session_types_2d =          ["of", "of_1m", "of_novel", "of_2m_novel", "of_dark", "ww", "mmaze"]
        session_types_navigation =  ["of", "of_1m", "of_novel", "of_2m_novel", "of_dark", "ww", "mmaze", "lt"]
        session_types_all =         ["of", "of_1m", "of_novel", "of_2m_novel", "of_dark", "ww", "mmaze", "lt", "sleep"]

        arena_names_2d =  ["of", "of_1m", "ww", "mmaze"];
        arena_names_all = ["of", "of_1m", "ww", "mmaze", "lt"];

        dirnames = struct( ...
            "lmt", "lmt_fits", ...
            "smdl", "smdl_fits")

        coverage_model_class = "ModularSweepAgentSimulation"
        coverage_model_default_sigma = 0.3
        coverage_model_default_kappa = 5

        alternation_score_chance_level = 0.4034; % tested with alternationScore4 2024-08-27

    end

    methods

        function self = SweepsSettings()
            %             self.fileDateStr = self.datestr(); % can be manually altered
            pcname = getenv("COMPUTERNAME");
            if isempty(pcname)
                %This will work on Linux
                pcname = getenv("HOSTNAME");
            end

            codeRoot = fileparts(fileparts(mfilename("fullpath")));
            self.codeRoot_ = codeRoot;

        end

        function val = cast(self, val)
            val = cast(val, self.numericClass);
            if self.useGpu
                val = gpuArray(val);
            end
        end

        function [a, str, units] = parseAngles(self, a, str)
            % Parse angles in *radians* into the required display format
            if nargin < 3 || isempty(str), str = ""; end
            units = self.angleUnits;
            if strcmpi(units, "degrees")
                a = rad2deg(a);
            elseif strcmpi(units, "cycles")
                a = a/(2*pi);
            end
            [str, units] = self.angleString(str);
        end

        function [str, units] = angleString(self, str)
            units = self.angleUnits;
            if strcmpi(units, "degrees")
                units = "°";
            elseif strcmpi(units, "radians")
                units = "rad.";
            elseif strcmpi(units, "cycles")
                units = "cycles";
            end
            str = str + " " + sprintf("(%s)", units);
        end

        function val = get.gve(self)
            val = self.gve_;

            % For all 2D pos ranges, generate "fine" and "coarse" binnings
            rngs = self.pos_lims;
            nbins = self.pos_nbins;
            rngFds = string(fieldnames(rngs))';
            binFds = string(fieldnames(nbins))';
            for rfd = rngFds
                for n = 1:2
                    if n==1
                        suffix = [];
                        coeff = 1;
                    else
                        suffix = "extended";
                        coeff = self.pos_extend_factor;
                    end
                    r = rngs.(rfd) * coeff;
                    for bfd = binFds
                        nb = nbins.(bfd) * coeff;
                        edges = linspace(r(1), r(2), nb+1)';
                        outFd = join(["pos", rfd, bfd, suffix], "_");
                        val.(outFd) = edges;
                    end
                end
            end


            for bfd = binFds
                nb = nbins.(bfd);
                edges = linspace(-pi, pi, nb+1)';
                outFd = "angular" + "_" + bfd;
                val.(outFd) = edges;
            end

            val = structfun(@(x) self.cast(x), val, "uni", 0);
        end

        function val = get.gv(self)
            val = structfun(@edg2cen, self.gve, "uni", 0);
        end

        function val = get.gverange(self)
            val = structfun(@(x) x([1, end]), self.gve, "uni", 0);
        end

        function val = get.gvspacing(self)
            val = structfun(@(x) x(2)-x(1), self.gv, "uni", 0);
        end

        function val = get.gg(self)
            fds = fieldnames(self.gve);
            for f = 1:numel(fds)
                fd = fds{f};
                bins = self.gv.(fd);
                if any(startsWith(fd, ["pos", "torus"]))
                    ndims = 2;
                else
                    ndims = 1;
                end
                if ndims == 1
                    val.(fd) = {bins};
                elseif ndims == 2
                    [xx, yy] = meshgrid(bins);
                    val.(fd) = {xx, yy};
                end
            end
        end

        function val = get.gg1(self)
            gg = self.gg;
            fds = fieldnames(self.gve);
            val = struct();
            for f = 1:numel(fds)
                fd = fds{f};
                vtmp = gg.(fd);
                vtmp = cellfun(@(x) x(:), vtmp, "uni", 0);
                val.(fd) = [vtmp{:}];
            end
        end

        function p = dataRoot(self)
            p = self.dataRoot_;
        end

        function p = codeRoot(self)
            p = self.codeRoot_;
            if isempty(p)
                userpath = getenv("USERPROFILE");
                p = fullfile(userpath, "OneDrive - NTNU", "Work", "2021_sweeps");
            end
        end

        function [p, parent] = filepath(self, varargin)
            % if self.useNewDataPath
                % p = fullfile(self.dataRootNew, varargin{:});
            % else
                p = fullfile(self.dataRoot, varargin{:});
            % end
            parent = fileparts(p);
        end

        function col = col_grid_module(self, n)
            col = [
                self.col_grid_m1
                self.col_grid_m2
                self.col_grid_m3
                self.col_grid_m4];
            if nargin == 2
                col = col(n, :);
            end
        end

        function cols = col_example_cell(self, n)
            cols = [
                0.85, 0, 0
                0, 0, 0
                1.0, 0.5, 0.0
                0.2, 0.2, 0.4
                0, 1, 0];
            if nargin == 2
                cols = cols(n, :);
            end
        end

        function nsamp = t2s(self, t, dt)
            % Time to samples
            if nargin < 3 || isempty(dt), dt = self.dt; end
            nsamp = t/dt;
        end

    end

end