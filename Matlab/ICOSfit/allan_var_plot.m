function [var_data,t_int]=allan_var_plot(data,time,varargin)

%function [var_data,t_int]=allan_var_plot(data,time)
%create allan variance plot given data and time)
%user is left to label axis, etc.
%Time is assumed to be in seconds.

max_t=(max(time)-min(time))/10;
min_t=.1;
max_tt=1;
inc=1;
j=1;
while max_t > max_tt
    min_t=min_t*10;
    if min_t==1;
        inc=1;
    elseif min_t==10;
        inc=5;
    elseif min_t==100;
        inc=10;
    elseif min_t==1000;
        inc=50;
    else
        inc=50;
    end
    max_tt=max_tt*10;
    if max_tt>max_t; max_tt=max_t; end
    for i=min_t:inc:max_tt
        [t_avg,data_avg]=binavg(time,data,i);
        var_data(j)=nanvar(data_avg);
        t_int(j)=i;
        j=j+1;
    end
end

%create plot
figure
nsubplot(3,1,[1,1],1);
plot(time-min(time),data)
axis tight
set(gca,'XAxisLocation','top')
xlabel('Time(s)')
box off
nsubplot(3,1,[3,2],1);
loglog(t_int,var_data)
hold on
x=[1,10,100,1e3,1e4];
var_i=nanvar(data);
%y=[var_i,var_i/sqrt(10),var_i/sqrt(100),var_i/sqrt(1e3),var_i/sqrt(1e4)];
y2=[var_i,var_i/10,var_i/100,var_i/1e3,var_i/1e4];
xlabel('Integration Time (s)')
ylabel('Allan Variance (\sigma^2)')
%plot(x,y,'--')
plot(x,y2,'-.')
if nargin==2
  text(500,var_i,['\sigma(1-sec) = ' num2str(nanstd(data))]);
  text(500,var_i/2,['\sigma(10-sec) = ' num2str(nanstd(fastavg(data,10)))]);
  text(500,var_i/4,['\sigma(100-sec) = ' num2str(nanstd(fastavg(data,100)))]);
elseif nargin==3
  text(500,var_i,['\sigma(1-sec) = ' num2str(nanstd(data))]);
  text(500,var_i/2,['\tau_{Allan} = ' num2str(varargin{1}) ' s']);
  text(500,var_i/4,['\sigma_{Allan} = ', num2str(nanstd(fastavg(data,varargin{1})))]);
end
box off