function showSegments(segment, averagingTimeInterval, divider)

fprintf('\n%s\nIndex\t\tStarting hour\t\tEnding hour\t\tAveraging interval: %d s (%.1f hr)\n%s', divider, averagingTimeInterval, averagingTimeInterval./3600, divider)

for idx = 1 : size(segment.stime, 2)
	fprintf('%d\t\t%10.1f\t\t%10.1f\n', idx, segment.stime(idx)./3600, segment.etime(idx)./3600)
end

fprintf('%s\n\n', divider)

end