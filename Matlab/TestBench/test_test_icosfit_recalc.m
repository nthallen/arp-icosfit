function aok = test_test_icosfit_recalc(verbose)
if nargin < 1; verbose = 0; end
warning off all;
warning off backtrace;
Summary = struct('success', 0, 'failure', 0, 'verbose', verbose);
Summary = bench(Summary, 'Basic1', {'ICOSout.R1.5p'}, '', ...
    'ICOSfit recalculation requires ICOSfit_format_ver >= 2');
Summary = bench(Summary, 'Basic2', {'ICOSout.R1.5p.R2.16rt', 'scans', 1:100});
Summary = bench(Summary, 'Basic3', {'ICOSout.R1.5p.R2.16rt', 1:100}, '', ...
    'Expected option string at argument 2');
Summary = bench(Summary, 'Basic4', {'ICOSout.R1.5p.R2.16rt', 'squirm'}, ...
    '', 'Unrecognized option: "squirm"');
Summary = bench(Summary, 'Basic5', {'ICOSout.R1.5p.R2.16rt', 'scans', 1:100});
Summary = bench(Summary, 'Basic6', ...
    {'ICOSout.R1.5p.R2.16rt', 'scans', 950:1005}, ...
    'Some scans in specified range are missing');
fprintf(1, '%d tests\n%d successes\n%d failures\n', ...
    Summary.success + Summary.failure, Summary.success, Summary.failure);
aok = Summary.failure == 0;

function Summary = bench(Summary, name, args, warnings, errors)
if nargin < 4
    warnings = '';
end
if nargin < 5
    errors = '';
end
lastwarn('');
success = 1;
msg = '';
try
    test_icosfit_recalc(args{:});
    [lw,li] = lastwarn;
    if ~strcmp(lw, warnings)
        success = 0;
        msg = sprintf( ...
            '%s: Expected warning "%s", received "%s"', ...
            name, warnings, lw);
    end
    if ~isempty(errors)
        success = 0;
        msg = sprintf( ...
            '%s: Expected error "%s", received none', ...
            name, errors);
    end
catch ME
    if ~strcmp(ME.message, errors)
        success = 0;
        msg = sprintf( ...
            '%s: Expected error "%s", received "%s"', ...
            name, errors, ME.message);
    end
end
if success
    Summary.success = Summary.success + 1;
else
    Summary.failure = Summary.failure + 1;
end
if Summary.verbose
    if success
        fprintf(1, '  OK: %s\n', name);
    else
        fprintf(1, 'FAIL: %s\n', msg);
    end
end
