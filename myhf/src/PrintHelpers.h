#pragma once


template <>
struct std::formatter<Eigen::MatrixXd>
{
	constexpr auto parse(std::format_parse_context& ctx) noexcept
	{
		return ctx.begin();
	}
	auto format(const Eigen::MatrixXd& mat, std::format_context& ctx) const
	{
		std::ostringstream oss;
		oss << std::setprecision(7) << mat;
		return std::format_to(ctx.out(), "{}", oss.str());
	}
};